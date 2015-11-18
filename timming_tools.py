#!/usr/bin/python
#####################################################
#   example.py - a program to ....                  #
#                                                   #
# Author: Dave Wheeler                              #
#                                                   #
# Purpose:                                          #
#                                                   #
# Usage:                                            #
#####################################################

from datetime import datetime
#would be nice to add some methods for seconds or minutes

def start_time():
    '''get a start time varible that can be used to time events'''
    return datetime.now()

def get_time(start_time):
    '''return the time since a start time either datetime.now or start_time()'''
    return datetime.now() - start_time

