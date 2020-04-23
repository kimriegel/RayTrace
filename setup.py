from distutils.core import setup, Extension

module1 = Extension('helloworld',
                    sources = ['demo.c']) # include_dirs=['/Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m/'])

setup (name = 'helloworld',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1])

