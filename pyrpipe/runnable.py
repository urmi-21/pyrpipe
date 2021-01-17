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
    
    Attributes
    ----------
    _runnable : bool
        Boolean indicating the runnable command
    _command : str
        Name of the unix command
    _param_yaml : str
        the name/path of .yaml file contating tool/command options (relative to the --param-dir)
    _args_style : str
        Style of commands: --key value(LINUX) or key=value(JAVA)
    _valid_args : dict or list
        A dict containing valid arguments for command subcommand, accessible as _valid_args[subcommand], or a list of valid options for command.
        
    """
    def __init__(self,*args,command=None,yaml=None,style='LINUX',deps=None,valid_args=None,**kwargs):
        """
        The runnable class to import any Unix command to python.

        Parameters
        ----------
        *args : tuple
            The positional arguments for the command. This will be saved as self._args
        command : String, optional
            The command name. The default is None.
        yaml : Str, optional
            The yaml file name containing command options. This file will be searched under the --param-dir value
        style : Str, optional
            The type of options passed --key value (LINUX) or key=value (JAVA). The default is 'LINUX'.
        deps : Str or List, optional
            A list of dependencies. These commands will be checked and exception is thrown if dependencies are not found. The default is None.
        valid_args : Dict or list, optional
            A dict containing the valid arguments of the command. This will be used to ignore invalid command options. The default is None.
        **kwargs : TYPE
            Command options passed to the tools. This will be saved as self._kwargs. These will override any commands supplied via the yaml file.

        Returns
        -------
        None.

        """
        self._runnable=True
        self._command=command
        self._deps=deps
        self._param_yaml=yaml
        #if not self._param_yaml and self._command:
        #    self._param_yaml=self._command+'.yaml'
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
        Str
            Return self.__param_yaml.

        """
        return self.__param_yaml
    @_param_yaml.setter
    def _param_yaml(self,value):
        """
        _param_yaml setter

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
        List
            Return self.__deps.

        """
        return self.__deps
    @_deps.setter
    def _deps(self,dep_list):
        """
        self._deps setter

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
        self._command setter

        Parameters
        ----------
        value : Str
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
        Resolve a tool parameter by passing as an argument.
        For example if unix command orfipy take a parameter --procs <num_threads>.
        This can be converted to a python variable accessible as an attribute of Runnable class.
        resolve_parameter function will update self._kwargs if parameter_keys exists. Otherwise it will create the parameter_key in self._kwargs
        To do this call the function as:
        <Runnable obj>.resolve_parameter("--procs",<passed_value>,<default_value>,'_threads')
        Now the <Runnable obj>._threads will point to "--procs" value

        Parameters
        ----------
        parameter_key : Str
            The parameter/option name for the Unix command/tool e.g. --threads
        passed_value : Str
            The value supplied by user.
        default_value : Str
            Default value to use if no value is supplied
        parameter_variable : Str
            Name of a variable that will be stored in the Runnable class. e.g. threads

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
        Check depndencies of a tool/command.

        Parameters
        ----------
        deps_list : List
            List of command to check.

        Raises
        ------
        OSError
            If a command is not found raise OSError.

        Returns
        -------
        bool
            Returns true is all commands are found.

        """
        if deps_list and not pe.check_dependencies(deps_list):
            #pu.print_boldred("ERROR. Please check dependencies for {}. Deps: {}".format(self._command," ".join(deps_list)))
            raise OSError("CommandNotFoundException")
        return True
                
    def load_args(self,*args,**kwargs):
        """
        Initializes the args (positonial arguments) and kwargs (options) passed during object creation.
        These are sored in self._args and self._kwargs for all future references.

        Parameters
        ----------
        *args : tuple
            Positional arguments
        **kwargs : dict
            the keyword arguments.

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
        Loads a .yaml file containing tool options/parameters

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
        Verify target file is present and is not LOCKED i.e. .Lock file is not present.

        Parameters
        ----------
        target : Str
            The target file.
        verbose : bool, optional
            Print additional messages. The default is False.

        Returns
        -------
        bool
            Return True is target is present and not Locked.

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
        Returns .Lock files associated with a target

        Parameters
        ----------
        target : Str
            Target file name.

        Returns
        -------
        lock_files : List
            List of .Lock files present.

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
        Verify a single target file is present

        Parameters
        ----------
        target : Str
            target file name.
        verbose : bool, optional
            Print additional messages. The default is False.

        Returns
        -------
        bool
            True is target file is presetn.

        """
        if not pu.check_files_exist(target):
            if verbose: pu.print_boldred("Error: target {} not found".format(target))
            return False
        return True
        
        
    
    def verify_target_list(self,target_list,verbose=False):
        """
        Verify a list of target files are present.

        Parameters
        ----------
        target_list : List
            List of target files.
        verbose : bool, optional
            Print additional messages. The default is False.

        Returns
        -------
        bool
            True is all targets are present.

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
        Cretes a temporary .Lock file associated with a target file and write a message in it.

        Parameters
        ----------
        target_list : List
            List of target files.
        message : Str
            Message to write in file.

        Returns
        -------
        templist : List
            A list of .Lock file names coressponding to the target files.

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
        Take a list of file names and removes them using os.remove

        Parameters
        ----------
        file_list : List
            List of file names ending with .Lock.

        Returns
        -------
        None.

        """
        for f in file_list:
            os.remove(f)
            
            
    
    def get_valid_parameters(self,subcommand):
        """
        Returns th evalid parameter list for a command subcommand after looking up the self_valid_args dictionary.

        Parameters
        ----------
        subcommand : Str
            The subcommand name

        Returns
        -------
        List
            Returns the list of valid options for the subcommand.

        """
        
        if isinstance(self._valid_args, list): return self._valid_args
        #if dict
        if isinstance(self._valid_args, dict):
            if subcommand and subcommand in self._valid_args:
                return self._valid_args[subcommand]
            
        #no valid format provided
        return None
        
    def run(self,*args, subcommand=None, target=None, requires=None, objectid=None, verbose=None,logs=None,**kwargs):
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
        #get valid args
        valid_args_subcommand=self.get_valid_parameters(subcommand)
        
        if subcommand:
            if isinstance(subcommand, str):
                subcommand=[subcommand]
            #add to command
            cmd.extend(subcommand)
            
        
        #parse and add parameters
        if self._args_style=='LINUX':
            cmd.extend(pu.parse_unix_args(valid_args_subcommand,kwargs))      
        elif self._args_style=='JAVA':
            cmd.extend(pu.parse_java_args(valid_args_subcommand,kwargs))
        else: 
            pu.print_boldred("Unknown args style: {}".format(self._args_style))
            raise ValueError("Unknown args style")
            
        
        #create locks on target; locks indicate incomplete commands
        if not _dryrun: locks=self.create_lock(target_list,' '.join(cmd))
            
        #execute command
        cmd_status=pe.execute_command(cmd,objectid=objectid,verbose=verbose,logs=logs)
        
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