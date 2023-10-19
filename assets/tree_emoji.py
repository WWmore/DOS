
"""https://github.com/earnestt1234/seedir
pip install seedir
pip install emoji
"""

import seedir as sd


#import os

# def foo(x): # omits folders with more than 2 items
#      if os.path.isdir(x) and len(os.listdir(x)) > 2:
#          return False
#      return True

folder = ['.git', '_pycache_', '__pycache__', '.spyproject', '.gitignore']
file = '.gitignore.txt'
path = r'C:\Users\WANGH0M\Documents\GitHub\DOS'
sd.seedir(path, depthlimit=2, exclude_folders=folder, exclude_files=file, style='emoji')