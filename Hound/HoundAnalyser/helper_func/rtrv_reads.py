import os


def rtrv_reads(path: str = None) -> list:
     """
          List files containing NGS illumina reads,
          and return a list with their absolute path.
     """
     if path is None:
         raise ValueError('I need path where illumina reads are stored.\n')
         os.system.exit(-1)
     else:
         dir_list = sorted(os.listdir(path))
         #dir_list.remove(".DS_Store")  # SH*T happens when coding in a mac...
         # Return only files, not directories
         rds_list = [os.path.abspath(path+direct) for direct in dir_list
                                                    if direct.find(".") > -1]
         return rds_list
