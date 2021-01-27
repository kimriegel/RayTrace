from distutils.core import setup, Extension

module1 = Extension('helloworld',
                    sources=['demo.c'])

setup(name='helloworld',
      version='1.0',
      description='This is a demo package',
      ext_modules=[module1])
