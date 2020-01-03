#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 18:20:41 2020

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os

class Quant:
    """This is an abstract class for quantification programs.
    """
    def __init__(self,index=""):
        self.category="Quantification"
        self.passedArgumentDict={}
        self.index=index
        
    def build_index(self):
        """function to create an index used by the quantification program
        """
        pass
    
    def check_index(self):
        """Function to check if index of this object is valid and exists
        """
    
    def perform_quant(self,sra_object):
        """Function to perform quant taking and sraobject as input
        
        """
        pass
    
