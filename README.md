FlybaseUpdater
==============

FlybaseUpdater is a Python Yapdi daemon to update modENCODE's data every so often with new gene information. It depends on [pyparallelcurl](https://github.com/petewarden/pyparallelcurl) and [YapDi](https://github.com/kasun/YapDi) to function.

To get started, simply edit `credentials.txt.example` to with your database's credentials and rename it to `credentials.txt`. After that, you can start FlybaseUpdater by doing `python FlybaseUpdater.py start` (though you might need sudo). You can also provide `restart` or `stop` in place of `start`.