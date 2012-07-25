# !/usr/bin/env python

''' FlybaseUpdater - Updates Flybase information for genes.
    Author - Shal Dengeki <shaldengeki@gmail.com>
    USAGE - python FlybaseUpdater.py start|stop|restart
	  REQUIRES - yapdi, MySQLdb, pyparallelcurl

    python FlybaseUpdater.py start would execute updateFlybase() in daemon mode 
    if there is no instance already running. 

    updateFlybase() connects to the MySQL database, loads a list of all the genes,
    and then hits all of the genes' Flybase pages, updating information as is necessary.

    python FlybaseUpdater.py stop would kill any running instance.

    python FlybaseUpdater.py restart would kill any running instance and
    start an instance. '''

import datetime
import httplib
import inspect
import os
import re
import subprocess
import sys
import syslog
import time
import traceback
import urllib
import urllib2

import yapdi
import MySQLdb
import pyparallelcurl

LOOP_EVERY_MINUTES = 15
PID_FILE = '/var/run/FlybaseUpdater.pid'
ISOFORM_IMAGE_ROOT_PATH = '/home/guoc/modENCODE/app/public/images/isoforms'
COMMAND_START = 'start'
COMMAND_STOP = 'stop'
COMMAND_RESTART = 'restart'

# causes httplib to return the partial response from a server in case the read fails to be complete.
def patch_http_response_read(func):
    def inner(*args):
        try:
            return func(*args)
        except httplib.IncompleteRead, e:
            return e.partial
    
    return inner
httplib.HTTPResponse.read = patch_http_response_read(httplib.HTTPResponse.read)

class DbConn(object):
  '''
  Simple database connection class to reconnect to MySQL if the connection times out.
  '''
  def __init__(self, username, password, database):
    self.username = username
    self.password = password
    self.database = database
    self.conn = None
    self.connect()
    self.cursor = self.conn.cursor()
  def connect(self):
    try:
      self.conn = MySQLdb.connect('localhost', self.username, self.password, self.database, charset="utf8", use_unicode=True)
    except MySQLdb.Error, e:
      print "Error connecting to MySQL database %d: %s" % (e.args[0],e.args[1])
      return False
    return True
  def queryDB(self, query, params=[], newCursor=False):
    try:
      if newCursor:
        cursor = self.conn.cursor()
      else:
        cursor = self.cursor
      cursor.execute(query, params)
    except (AttributeError, MySQLdb.OperationalError):
      # lost connection. reconnect and re-query.
      if not self.connect():
        print "Unable to reconnect to MySQL."
        return False
      cursor = self.conn.cursor()
      cursor.execute(query, params)
      self.cursor = cursor
    self.commit()
    return cursor
  def commit(self):
    self.conn.commit()

def usage():
  print("USAGE: python %s %s|%s|%s" % (sys.argv[0], COMMAND_START, COMMAND_STOP, COMMAND_RESTART))

def downloadGeneIsoformImage(image_url, dest_url, gene_id, dbConn):
  try:
    isoform_image = urllib2.urlopen(image_url)
  except urllib2.HTTPError:
    pass
  else:
    CHUNK = 16 * 1024
    with open(dest_url, 'wb') as fp:
      while True:
        chunk = isoform_image.read(CHUNK)
        if not chunk: break
        fp.write(chunk)
      foo = dbConn.queryDB(u'''UPDATE `transcription_factors` SET `isoform_image_path` = %s WHERE `id` = %s''', [str(os.path.basename(dest_url)), str(gene_id)], newCursor=True)

  
def compareGeneFlybaseInfo(flybase_html, flybase_url, ch, params):
  geneDBInfo = params['geneDBInfo']
  dbConn = params['dbConn']

  # construct a hash of aliases consisting of the common name, gene symbol name, and any "also known as" strings specified on the page.
  geneInfo = {}
  geneInfo['aliases'] = {}
  symbol = re.search('Symbol\<\/th\>\<td\>\<span\ class\=\"greytext\"\>Dmel\\\<\/span\>(.*?)\<\/td\>', flybase_html)
  if symbol and len(symbol.groups()) > 0:
    geneInfo['name'] = symbol.groups()[0].strip()
    geneInfo['aliases'][geneInfo['name']] = 1
  akas = re.search('Also\ Known\ As\<\/th\>\<td\ colspan\=\"3\"\>(.*?)\<\/td\>', flybase_html)
  if akas and len(akas.groups()) > 0:
    for item in akas.groups()[0].split(", "):
      geneInfo['aliases'][item.strip()] = 1
  # compare these aliases against the ones in the DB.
  aliasesToInsert = []
  aliasesToDelete = []
  for aliasName in geneInfo['aliases'].keys():
    if aliasName not in geneDBInfo['aliases']:
      aliasesToInsert.extend([str(aliasName), str(geneDBInfo['id'])])
  for aliasName in geneDBInfo['aliases'].keys():
    if aliasName not in geneInfo['aliases']:
      aliasesToDelete.extend(str(geneDBInfo['aliases'][aliasName]))
  aliasInsertQuery = []
  for x in range(len(aliasesToInsert)/2):
    aliasInsertQuery.append("(%s, %s)")
  if len(aliasInsertQuery) > 0:
    foo = dbConn.queryDB(u'''INSERT IGNORE INTO `aliases` (`name`, `transcription_factor_id`) VALUES ''' + ",".join(aliasInsertQuery), aliasesToInsert, newCursor=True)
    # syslog.syslog(syslog.LOG_NOTICE, "Inserted " + str(len(aliasesToInsert)/2) + " aliases.")
  
  if len(aliasesToDelete) > 0:
    foo = dbConn.queryDB(u'''DELETE FROM `aliases` WHERE `id` IN (''' + ",".join(['%s' for x in aliasesToDelete]) + ''')''', aliasesToDelete, newCursor=True)
    # syslog.syslog(syslog.LOG_NOTICE, "Deleted " + str(len(aliasesToDelete)) + " aliases.")
  
  # construct a hash of hashes of isoforms and their external IDs.
  geneInfo['isoforms'] = {}
  gene_isoform_table = re.search('\(aa\)</div></div></div></div></div>(.*?)</td></tr>', flybase_html, flags=re.DOTALL)
  if gene_isoform_table and len(gene_isoform_table.groups()) > 0:
    for item in gene_isoform_table.groups()[0].split("<div class=\"line-wrapper\">"):
      isoform_name = re.search('<div\ class\=\"trans\_name\"><a\ href=\"\/reports\/[0-9a-zA-Z]+\.html\">(<span\ style\=\"background\-color\:\#bbffbb\;\">)?(?P<name>.*?)(</span>)?</a></div>', item)
      if isoform_name:
        isoform = {'name': isoform_name.group("name").strip(), 'flybase_id': '', 'refseq_id': ''}
        isoform_flybase = re.search('<div\ class\=\"trans\_first\_block\"><div\ class\=\"trans\_ID\">(?P<flybase_id>[0-9A-Za-z]+)</div><div\ class\=\"trans\_second\_block\">', item)
        if isoform_flybase:
          isoform['flybase_id'] = isoform_flybase.group("flybase_id").strip()
        isoform_refseq = re.search('<a\ href\=\"http\:\/\/www\.ncbi\.nlm\.nih\.gov\/entrez\/viewer\.fcgi\?tool\=FlyBase\&amp\;val\=(?P<refseq_id>[0-9A-Za-z_]+)\">', item)
        if isoform_refseq:
          isoform['refseq_id'] = isoform_refseq.group("refseq_id").strip()
        geneInfo['isoforms'][isoform['name']] = isoform
  # compare these isoforms against the ones in the DB.
  isoformsToInsert = []
  isoformsToDelete = []
  for isoformName in geneInfo['isoforms'].keys():
    if isoformName not in geneDBInfo['isoforms']:
      isoformsToInsert.extend([str(isoformName), str(geneDBInfo['id']), str(geneInfo['isoforms'][isoformName]['flybase_id']), str(geneInfo['isoforms'][isoformName]['refseq_id'])])
  for isoformName in geneDBInfo['isoforms'].keys():
    if isoformName not in geneInfo['isoforms']:
      isoformsToDelete.extend(str(geneDBInfo['isoforms'][isoformName]['id']))
    else:
      geneInfo['isoforms'][isoformName]['id'] = geneDBInfo['isoforms'][isoformName]['id']
      if geneInfo['isoforms'][isoformName] != geneDBInfo['isoforms'][isoformName]:
        # update the DB with Flybase's info.
        foo = dbConn.queryDB(u'''UPDATE `isoforms` SET `flybase_id` = %s, `refseq_id` = %s WHERE `id` = %s LIMIT 1''', [str(geneInfo['isoforms'][isoformName]['flybase_id']), str(geneInfo['isoforms'][isoformName]['refseq_id']), str(geneDBInfo['isoforms'][isoformName]['id'])], newCursor=True)
  isoformInsertQuery = []
  for x in range(len(isoformsToInsert)/4):
    isoformInsertQuery.append("(%s, %s, %s, %s)")
  if len(isoformInsertQuery) > 0:
    foo = dbConn.queryDB(u'''INSERT IGNORE INTO `isoforms` (`name`, `transcription_factor_id`, `flybase_id`, `refseq_id`) VALUES ''' + ",".join(isoformInsertQuery), isoformsToInsert, newCursor=True)
    # syslog.syslog(syslog.LOG_NOTICE, "Inserted " + str(len(isoformsToInsert)/4) + " isoforms.")
  
  if len(isoformsToDelete) > 0:
    foo = dbConn.queryDB(u'''DELETE FROM `isoforms` WHERE `id` IN (''' + ",".join(['%s' for x in isoformsToDelete]) + ''')''', isoformsToDelete, newCursor=True)
    # syslog.syslog(syslog.LOG_NOTICE, "Deleted " + str(len(isoformsToDelete)) + " isoforms.")
  
  # get isoform image size and url.
  geneInfo['isoform_image'] = {'url': '', 'size': 0}
  geneDBInfo['isoform_image'] = {'url': os.path.join(ISOFORM_IMAGE_ROOT_PATH, geneDBInfo['isoform_image_path']), 'size': 0}
  isoform_image_path = re.search('class\=\"noborder\-line\-wrapper\"><img\ align\=\"middle\"\ alt\=\"detailed view\"\ width\=\"[0-9]+\"\ name\=\"detailedView\"\ height\=\"[0-9]+\"\ border\=\"0\"\ src\=\"(.*?)\"\ usemap', flybase_html)
  if isoform_image_path:
    geneInfo['isoform_image']['url'] = "".join(['http://flybase.org',isoform_image_path.groups()[0]])
    isoform_image = urllib.urlopen(geneInfo['isoform_image']['url'])
    isoform_image_meta = isoform_image.info()
    geneInfo['isoform_image']['size'] = isoform_image_meta.getheaders("Content-Length")[0]
    # compare the image filenames and if necessary re-download.
    if os.path.basename(geneInfo['isoform_image']['url']) != geneDBInfo['isoform_image_path']:
      try:
        os.remove(geneDBInfo['isoform_image']['url'])
      except (IOError, OSError):
        pass
      geneDBInfo['isoform_image']['url'] = os.path.join(ISOFORM_IMAGE_ROOT_PATH, os.path.basename(geneInfo['isoform_image']['url']))
      downloadGeneIsoformImage(geneInfo['isoform_image']['url'], geneDBInfo['isoform_image']['url'], geneDBInfo['id'], dbConn)
      # syslog.syslog(syslog.LOG_NOTICE, "Downloaded an updated isoform image.")
    else:
      # get isoform image size on disk (if one exists).
      try:
        f = open(geneDBInfo['isoform_image']['url'], "rb")
      except (IOError, OSError):
        geneDBInfo['isoform_image']['size'] = 0
      else:
        geneDBInfo['isoform_image']['size'] = len(f.read())
        f.close()
      if geneDBInfo['isoform_image']['size'] != geneInfo['isoform_image']['size']:
        try:
          os.remove(geneDBInfo['isoform_image']['url'])
        except (IOError, OSError):
          pass
        geneDBInfo['isoform_image']['url'] = os.path.join(ISOFORM_IMAGE_ROOT_PATH, os.path.basename(geneInfo['isoform_image']['url']))
        downloadGeneIsoformImage(geneInfo['isoform_image']['url'], geneDBInfo['isoform_image']['url'], geneDBInfo['id'], dbConn)
        # syslog.syslog(syslog.LOG_NOTICE, "Downloaded an updated isoform image.")


def updateFlybase():
  '''
  Updates genes in a MySQL database with Flybase info.
  '''
  syslog.openlog("FlybaseUpdater.info", 0, syslog.LOG_USER)
  try:
    dbConn = DbConn(MYSQL_USERNAME, MYSQL_PASSWORD, MYSQL_DATABASE)
    if not dbConn:
      syslog.syslog(syslog.LOG_NOTICE, "Unable to connect to MySQL.")
      return
      
    # Now loop every N minutes, updating flybase info.
    while 1:
      try:
        loopStartTime = time.time()
        parallelcurl = pyparallelcurl.ParallelCurl(20)
        # load all genes (along with alias lists) into memory.
        geneCursor = dbConn.queryDB(u'''SELECT `id`, `name`, `flybase_id`, `cg_id`, `isoform_image_path` FROM `transcription_factors` WHERE (`flybase_id` IS NOT NULL && `flybase_id` != '')''')
        gene = geneCursor.fetchone()
        while gene is not None:
          # assemble a list of aliases for this gene.
          geneInfo = {'id': int(gene[0]), 'name': str(gene[1]), 'flybase_id': str(gene[2]), 'cg_id': str(gene[3]), 'isoform_image_path': str(gene[4]), 'aliases': {}, 'isoforms': {}}
          geneAliasQuery = dbConn.queryDB(u'''SELECT `id`, `name` FROM `aliases` WHERE `transcription_factor_id` = %s''', params=[str(gene[0])], newCursor=True)
          alias = geneAliasQuery.fetchone()
          while alias is not None:
            geneInfo['aliases'][str(alias[1])] = int(alias[0])
            alias = geneAliasQuery.fetchone()
          
          # assemble a list of isoforms for this gene.
          geneIsoformQuery = dbConn.queryDB(u'''SELECT `id`, `name`, `flybase_id`, `refseq_id` FROM `isoforms` WHERE `transcription_factor_id` = %s''', params=[str(gene[0])], newCursor=True)
          isoform = geneIsoformQuery.fetchone()
          while isoform is not None:
            geneInfo['isoforms'][isoform[1].encode('utf-8')] = {'id': int(isoform[0]), 'name': isoform[1].encode('utf-8'), 'flybase_id': isoform[2].encode('utf-8'), 'refseq_id': isoform[3].encode('utf-8')}
            isoform = geneIsoformQuery.fetchone()
          
          # submit this request to parallelcurl object to be compared.
          parallelcurl.startrequest('http://flybase.org/reports/'+ str(gene[2]) +'.html', compareGeneFlybaseInfo, {'geneDBInfo': geneInfo, 'dbConn': dbConn})
          gene = geneCursor.fetchone()
        parallelcurl.finishallrequests()
        syslog.syslog(syslog.LOG_NOTICE, "Loop finished. Sleeping for " + str(LOOP_EVERY_MINUTES) + " minutes.")
      except:
        syslog.syslog(syslog.LOG_NOTICE, "Recoverable error:\n" + str(traceback.format_exc()) + "\n")
      time.sleep(LOOP_EVERY_MINUTES * 60)
  except:
    syslog.syslog(syslog.LOG_NOTICE, "Fatal error:\n" + str(traceback.format_exc()) + "\n")

# Invalid executions
if len(sys.argv) < 2 or sys.argv[1] not in [COMMAND_START, COMMAND_STOP, COMMAND_RESTART]:
  usage()
  exit()

# Load credentials from a textfile.
openCredentialsFile = open('./credentials.txt')
mysqlLogin = openCredentialsFile.readline().strip('\n').split(',')
if len(mysqlLogin) < 2:
  print "MySQL login not found in credentials file."
  exit()

MYSQL_USERNAME = mysqlLogin[0].strip()
MYSQL_PASSWORD = mysqlLogin[1].strip()
MYSQL_DATABASE = openCredentialsFile.readline().strip()
if len(MYSQL_DATABASE) < 1:
  print "Database data not found in credentials file."
  exit()
openCredentialsFile.close()

if sys.argv[1] == COMMAND_START:
  daemon = yapdi.Daemon(pidfile=PID_FILE)

  # Check whether an instance is already running
  if daemon.status():
    print("An instance is already running.")
    exit()
  retcode = daemon.daemonize()
  # Execute if daemonization was successful else exit
  if retcode == yapdi.OPERATION_SUCCESSFUL:
    updateFlybase()
  else:
    syslog.syslog(syslog.LOG_NOTICE, 'Daemonization failed')

elif sys.argv[1] == COMMAND_STOP:
  daemon = yapdi.Daemon(pidfile=PID_FILE)

  # Check whether no instance is running
  if not daemon.status():
    print("No instance running.")
    exit()
  retcode = daemon.kill()
  if retcode == yapdi.OPERATION_FAILED:
    print('Trying to stop running instance failed')

elif sys.argv[1] == COMMAND_RESTART:
  daemon = yapdi.Daemon(pidfile=PID_FILE)
  retcode = daemon.restart()
  # Execute if daemonization was successful else exit
  if retcode == yapdi.OPERATION_SUCCESSFUL:
    updateFlybase()
  else:
    print('Daemonization failed')