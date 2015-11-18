"""
log.py
=====

Setup logging utils for nested module logging

Adapted from the accepted answer here: http://stackoverflow.com/questions/7621897/python-logging-module-globally
"""

import logging

def createCustomLogger(name):
    formatter = logging.Formatter(fmt='[%(asctime)s][%(levelname)s][%(module)s] %(message)s', datefmt='%m/%d %I:%M:%S%p')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger
