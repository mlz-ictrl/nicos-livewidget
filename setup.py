from distutils.core import setup, Extension
import os
from os import path
import platform
from sipdistutils import build_ext as sip_build_ext
from PyQt4 import pyqtconfig

config = pyqtconfig.Configuration()
pyqt_sip_flags = config.pyqt_sip_flags.split()
pyqwt_sip_flags = ['-t', 'Qwt_5_2_0']


class moc_build_ext(sip_build_ext):
    '''Build with moc-generated files '''

    def finalize_options(self):
        sip_build_ext.finalize_options(self)
        if isinstance(self.sip_opts, str):
            self.sip_opts = self.sip_opts.split(' ')

        self.sip_opts = self.sip_opts + pyqt_sip_flags + pyqwt_sip_flags

    def _sip_sipfiles_dir(self):
        return pyqtconfig.Configuration().pyqt_sip_dir

    def swig_sources(self, sources, extension=None):
        # Create .moc files from headers
        ret = sip_build_ext.swig_sources(self, sources, extension)
        moc_bin = config.build_macros()['MOC']
        for source in sources:
            if not source.endswith('.cpp'):
                continue
            header = source[:-4] + '.h'
            if not path.isfile(header):
                continue
            dest = path.join(self.build_temp, 'moc_' + source)
            self.spawn([moc_bin] + ['-o', dest, header])
            if path.getsize(dest) == 0:
                continue
            ret.append(dest)
        return ret


dist = platform.linux_distribution()[0].strip()  # old openSUSE appended a space here :(
if dist == 'openSUSE':
    extra_include_dirs = ["/usr/include/qwt5",
                          "/usr/include/qwt",
                          "/usr/include/libcfitsio0",
                          "/usr/include/cfitsio"]
    extra_libs = ["qwt", "cfitsio"]
elif dist in ['Ubuntu', 'LinuxMint']:
    extra_include_dirs = ["/usr/include/qwt-qt4"]
    extra_libs = ["qwt-qt4", "cfitsio"]
elif dist == 'CentOS':
    extra_include_dirs = ["/usr/local/qwt5/include"]
    extra_libs = ["qwt", "cfitsio"]
else:
    print("WARNING: Don't know where to find Qwt headers and libraries "
          "for your distribution")
    # still try to build with usable defaults
    extra_include_dirs = ["/usr/include/qwt5", "/usr/include/qwt"]
    extra_libs = ["qwt", "cfitsio"]

extra_include_dirs.extend(path.join(config.qt_inc_dir, subdir)
                          for subdir in ['', 'QtCore', 'QtGui'])
extra_libs.extend(['QtCore', 'QtGui'])
extra_lib_dirs = [config.qt_lib_dir]

sources = [cppfile for cppfile in os.listdir('.')
           if cppfile.startswith('lw_') and cppfile.endswith('.cpp')
           and cppfile != 'lw_main.cpp']

setup(
    name='nicoslivewidget',
    version='1.0',
    ext_modules=[
        Extension('nicoslivewidget',
                  ['livewidget.sip'] + sources,
                  include_dirs=['.'] + extra_include_dirs,
                  library_dirs=extra_lib_dirs,
                  libraries=extra_libs,
                  ),
    ],
    cmdclass={'build_ext': moc_build_ext}
)
