#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 13:55:13 2020

@author: usingh
"""
import os
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
from pyrpipe import param_loader as pl
from pyrpipe import _params_dir
from pyrpipe import _dryrun
from pyrpipe import _force
import tempfile 

class Runnable:
    """The runnable class
    """
    def __init__(self,*args,command=None,yaml=None,style='LINUX',deps=None,valid_args=None,**kwargs):
        """
        

        Parameters
        ----------
        *args : TYPE
            DESCRIPTION.
        command : TYPE, optional
            DESCRIPTION. The default is None.
        yaml : TYPE, optional
            DESCRIPTION. The default is None.
        style : TYPE, optional
            DESCRIPTION. The default is 'LINUX'.
        deps : TYPE, optional
            DESCRIPTION. The default is None.
        valid_args : TYPE, optional
            DESCRIPTION. The default is None.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self._runnable=True
        self._command=command
        self._deps=deps
        self._param_yaml=yaml
        self._args_style=style
        #valid_args can be None or a list, or a dict if subcommands are used
        self._valid_args=valid_args
        self.load_args(*args,**kwargs)
        self.load_yaml()
            
    
    @property
    def _param_yaml(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.__param_yaml
    @_param_yaml.setter
    def _param_yaml(self,value):
        """
        

        Parameters
        ----------
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.__param_yaml=value
        self.load_yaml()
            
            
    @property
    def _deps(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.__deps
    @_deps.setter
    def _deps(self,dep_list):
        """
        

        Parameters
        ----------
        dep_list : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #validate command exists
        if dep_list and self.check_dependency(dep_list):
            self.__deps=dep_list
     
        
    @property
    def _command(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.__command
    @_command.setter
    def _command(self,value):
        """
        

        Parameters
        ----------
        value : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #if already exists
        if hasattr(self,'_command') and self._command:
            raise TypeError("Can not modify _command")
        
        #validate command exists
        if value and self.check_dependency([value]):
            self.__command=value        
        
    
    def resolve_parameter(self,parameter_key,passed_value,default_value,parameter_variable):
        """
        

        Parameters
        ----------
        parameter_key : TYPE
            DESCRIPTION.
        passed_value : TYPE
            DESCRIPTION.
        default_value : TYPE
            DESCRIPTION.
        parameter_variable : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #First check the passed value; then the passed **kwargs
        #None values will be ignored
        if passed_value:
            setattr(self, parameter_variable, str(passed_value))
            self._kwargs[parameter_key]=str(passed_value)
        elif parameter_key in self._kwargs and self._kwargs[parameter_key]:
            setattr(self, parameter_variable, self._kwargs[parameter_key])
        elif default_value:
            setattr(self, parameter_variable, str(default_value))
            self._kwargs[parameter_key]=str(default_value)
         
    def check_dependency(self,deps_list):
        """
        

        Parameters
        ----------
        deps_list : TYPE
            DESCRIPTION.

        Raises
        ------
        OSError
            DESCRIPTION.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        if deps_list and not pe.check_dependencies(deps_list):
            #pu.print_boldred("ERROR. Please check dependencies for {}. Deps: {}".format(self._command," ".join(deps_list)))
            raise OSError("CommandNotFoundException")
        return True
                
    def load_args(self,*args,**kwargs):
        """
        

        Parameters
        ----------
        *args : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #init the parameters for the object
        if args:
            self._args=args
        else:
            self._args=()
        if kwargs:
            self._kwargs=kwargs
        else:
            self._kwargs={}
        
    def load_yaml(self):        
        """
        

        Returns
        -------
        None.

        """
        #read yaml parameters
        if not self._param_yaml:
            return
        yamlfile=os.path.join(_params_dir,self._param_yaml)
        #override yaml parameters by kwargs
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            self._kwargs={**yaml_kwargs,**self._kwargs}
            
            
    def verify_integrity(self,target,verbose=False):
        """
        

        Parameters
        ----------
        target : TYPE
            DESCRIPTION.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        #check if lock exists
        lock_files=self.get_lock_files(target)
        if len(lock_files) >0:
            #remove the target and locks
            if pu.check_files_exist(target):
                pu.print_notification("Found incomplete target {}. Restarting command...".format(target))
                self.remove_locks(lock_files+[target])
            else: self.remove_locks(lock_files)
        
        return True
    
    def get_lock_files(self,target):
        """
        

        Parameters
        ----------
        target : TYPE
            DESCRIPTION.

        Returns
        -------
        lock_files : TYPE
            DESCRIPTION.

        """
        #check if lock exists
        filepath=pu.get_file_directory(target)
        filename=pu.get_filename(target)
        pre='.*'
        suff='.*\.Lock$'
        pattern=pre+filename+suff
        lock_files=pu.find_files(filepath,pattern)
        return lock_files
        
            
    def verify_target(self,target,verbose=False):
        """
        

        Parameters
        ----------
        target : TYPE
            DESCRIPTION.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        if not pu.check_files_exist(target):
            if verbose: pu.print_boldred("Error: target {} not found".format(target))
            return False
        return True
        
        
    
    def verify_target_list(self,target_list,verbose=False):
        """
        

        Parameters
        ----------
        target_list : TYPE
            DESCRIPTION.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        #nothing to verify
        if not target_list:
            return True
        for target in target_list:
            if not self.verify_target(target,verbose):
                return False
        return True
    
    
    def create_lock(self,target_list,message):
        """
        

        Parameters
        ----------
        target_list : TYPE
            DESCRIPTION.
        message : TYPE
            DESCRIPTION.

        Returns
        -------
        templist : TYPE
            DESCRIPTION.

        """
        templist=[]
        for f in target_list:
            temp_path=pu.get_file_directory(f)
            if not pu.check_paths_exist(temp_path): pu.mkdir(temp_path)
            prefix=pu.get_filename(f)+'_'
            temp = tempfile.NamedTemporaryFile(prefix=prefix,suffix='.Lock', dir=temp_path,delete=False)
            #TODO: dump command in lock
            timestamp=pu.get_timestamp()
            temp.write(str.encode(timestamp+'\t'+message))

            templist.append(temp.name)
        return templist
    
    def remove_locks(self,file_list):
        """
        

        Parameters
        ----------
        file_list : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for f in file_list:
            os.remove(f)
            
            
    
    def get_valid_parameters(self,subcommand):
        """
        

        Parameters
        ----------
        subcommand : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if subcommand:
            return self._valid_args[subcommand]
        return self._valid_args
        pass
        
    def run(self,*args, subcommand=None, target=None, requires=None, objectid=None, **kwargs):
        """
        

        Parameters
        ----------
        *args : Tuple
            Positoinal arguments passed to a command. This will copmletely REPLACE the exsiting self._args created during initialization of the runnable object.
        subcommand : String or List, optional
            DESCRIPTION. subcommand passed to the command. The default is None.
        target : Str or List of Str, optional
            DESCRIPTION. The expected output/target files produced by the run operation. False  is returned is all target files are not found after the command. The default is None.
        requires : Str or List of Str, optional
            DESCRIPTION. Files required to strat the run method. Exception is thrown if files are missing. The default is None.
        objectid : Str, optional
            DESCRIPTION. A uniq id to identify the run operation in the logs. Thi is useful for benchmarks. The default is None.
        **kwargs : Keyword arguments
            DESCRIPTION. The options to be passed to the command. This will OVERRIDE ANY EXISTING options in the self._kwargs created during initialization of the runnable object.

        Raises
        ------
        TypeError
            If incorerct types are used for target and required.
        FileNotFoundError
            Raises FileNotFoundError if any of the required files are missing.
        OSError
            Raises OSError if the command is incorrect or not present in path.
        ValueError
            Raises ValueError if args_type is something other than LINUX or JAVA.

        Returns
        -------
        bool
            Return the status of command as True or False. True implies command had 0 exit-code and all target files were found after the command finished.

        """
        
        #create target list
        target_list=[]
        locks=[]
        requires_list=[]
        
        if target:
            if isinstance(target, str):
                target_list=[target]
            elif isinstance(target, list):
                target_list=target
            else:
                raise TypeError("target must be a string or a list object")
        
        #ckeck for locks and remove previous locks and associated targets if exist
        for target in target_list:
            self.verify_integrity(target)        
        
        #if target already present and not overwrite exists then return
        if not _force and target_list:
            if self.verify_target_list(target_list):
                pu.print_green('Target files {} already exist.'.format(', '.join(target_list)))
                return True
            
        #check if all requirements are satisfied
        if requires:
            if isinstance(requires, str):
                requires_list=[requires]
            elif isinstance(requires, list):
                requires_list=requires
            else:
                raise TypeError("requires must be a string or a list object")
        
        #Raise exception if requirements not satisfied
        if requires_list:
            if not self.verify_target_list(requires_list):
                pu.print_boldred('Required files {} fot found.'.format(', '.join(requires_list)))
                raise FileNotFoundError("FilesNotFound")
            #check if any required file had lock
            for file in requires_list:
                if len(self.get_lock_files(file)):
                    pu.print_boldred('Required file {} is corrupt. Please verify file is correct and remove any .Lock files'.format(', '.join(requires_list)))
                    raise FileNotFoundError("FilesNotFound")        
        
        
        #override class kwargs by passed kwargs
        kwargs={**self._kwargs,**kwargs}
        #if no args provided use constructor's args
        if not args:
            args=self._args
        #if args are not None
        if args and args[0]:
            kwargs['--']=args
            
        
        #make a copy of self._command
        if not self._command:
            pu.print_boldred("Error: command can not be None or empty")
            raise OSError("CommandNotFoundException")    
            
        cmd=[]
        if isinstance(self._command, list):
            cmd=self._command.copy()            
        elif isinstance(self._command, str):
            cmd=[self._command]
        
        #if subcommand supplied
        if subcommand:
            if isinstance(subcommand, str):
                subcommand=[subcommand]
            #add to command
            cmd.extend(subcommand)
            
            
        #parse and add parameters
        if self._args_style=='LINUX':
            cmd.extend(pu.parse_unix_args(self._valid_args,kwargs))      
        elif self._args_style=='JAVA':
            cmd.extend(pu.parse_java_args(self._valid_args,kwargs))
        else: 
            pu.print_boldred("Unknown args style: {}".format(self._args_style))
            raise ValueError("Unknown args style")
            
        
        #create locks on target; locks indicate incomplete commands
        if not _dryrun: locks=self.create_lock(target_list,' '.join(cmd))
            
        #execute command
        cmd_status=pe.execute_command(cmd,objectid=objectid)
        
        # if command finished remove locks
        self.remove_locks(locks)
        
        if not cmd_status:
            pu.print_boldred("{} failed: {}".format(self._command," ".join(cmd)))
            #remove target files
            if not _dryrun and target_list:
                pu.print_boldred("Removing target files {}: ".format(', '.join(target_list)))
                pe.delete_files(*target_list)
            return False
     
        
        if cmd_status and target_list and not _dryrun:
            return self.verify_target_list(target_list,verbose=True)
        
        #return status
        return cmd_status