import os
import sys
import platform
from os import path
from setuptools import setup, Extension
from distutils.command.config import config as CommandConfig
from distutils.dist import Distribution

from sipdistutils import build_ext as sip_build_ext
from gitversion import get_git_version

try:
    from PyQt4 import pyqtconfig
except ImportError:
    from sipconfig import Configuration

    def query_var(varname):
        p = os.popen('qmake -query ' + varname, 'r')
        return p.read().strip()
    # new-style configured PyQt4, no pyqtconfig module
    from PyQt4.QtCore import PYQT_CONFIGURATION
    pyqt_sip_flags = PYQT_CONFIGURATION['sip_flags'].split()
    pyqt_sip_dir = path.join(Configuration().default_sip_dir, 'PyQt4')
    moc_bin = path.join(query_var('QT_INSTALL_BINS'), 'moc')
    qt_inc_dir = query_var('QT_INSTALL_HEADERS')
    qt_lib_dir = query_var('QT_INSTALL_LIBS')
else:
    config = pyqtconfig.Configuration()
    pyqt_sip_flags = config.pyqt_sip_flags.split()
    moc_bin = config.build_macros()['MOC']
    pyqt_sip_dir = config.pyqt_sip_dir
    qt_inc_dir = config.qt_inc_dir
    qt_lib_dir = config.qt_lib_dir

pyqwt_sip_flags = ['-t', 'Qwt_5_2_0']


class moc_build_ext(sip_build_ext):
    '''Build with moc-generated files '''

    def finalize_options(self):
        sip_build_ext.finalize_options(self)
        if isinstance(self.sip_opts, str):
            self.sip_opts = self.sip_opts.split(' ')

        self.sip_opts = self.sip_opts + pyqt_sip_flags + pyqwt_sip_flags
        get_git_version()  # create version.h and RELEASE-VERSION files

    def _sip_sipfiles_dir(self):
        return pyqt_sip_dir

    def swig_sources(self, sources, extension=None):
        # Create .moc files from headers
        ret = sip_build_ext.swig_sources(self, sources, extension)
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


if sys.platform == 'darwin':
    extra_include_dirs = ["/usr/local/qwt/include", "/usr/local/cfitsio/include"]
    extra_libs = ["qwt", "cfitsio", "tiff"]
    extra_lib_dirs = ["/usr/local/qwt/lib", "/usr/local/cfitsio/lib"]
else:
    extra_lib_dirs = []
    dist = platform.linux_distribution()[0].strip()  # old openSUSE appended a space here :(
    if dist == 'openSUSE':
        extra_include_dirs = ["/usr/include/qwt5",
                              "/usr/include/qwt",
                              "/usr/include/libcfitsio0",
                              "/usr/include/cfitsio"]
        extra_libs = ["qwt", "cfitsio", "tiff"]
    elif dist in ['Ubuntu', 'LinuxMint', 'debian',]:
        extra_include_dirs = ["/usr/include/qwt-qt4", "/usr/include/qwt"]
        extra_libs = ["qwt-qt4", "cfitsio", "tiff"]
    elif dist == 'CentOS':
        extra_include_dirs = ["/usr/local/qwt5/include"]
        extra_libs = ["qwt", "cfitsio", "tiff"]
        extra_lib_dirs.append("/usr/local/qwt5/lib")
    elif dist == 'Fedora':
        extra_include_dirs = ["/usr/include/qwt5-qt4", "/usr/include/cfitsio"]
        extra_libs = ["qwt5-qt4", "cfitsio", "tiff"]
    else:
        print("WARNING: Don't know where to find Qwt headers and libraries "
              "for your distribution")
        # still try to build with usable defaults
        extra_include_dirs = ["/usr/include/qwt5", "/usr/include/qwt"]
        extra_libs = ["qwt", "cfitsio", "tiff"]

extra_include_dirs.extend(path.join(qt_inc_dir, subdir)
                          for subdir in ['', 'QtCore', 'QtGui'])
extra_libs.extend(['QtCore', 'QtGui'])
extra_lib_dirs.append(qt_lib_dir)

sources = [cppfile for cppfile in os.listdir('.')
           if cppfile.startswith('lw_') and cppfile.endswith('.cpp')
           and cppfile != 'lw_main.cpp']

conf = CommandConfig(Distribution())

devfits = conf.check_header('fitsio.h', extra_include_dirs)
devtiff = conf.check_header('tiff.h', extra_include_dirs)
if not devfits or not devtiff:
    if not devfits:
        conf.warn("Please install developer files of the 'fitsio' library.")
    if not devtiff:
        conf.warn("Please install developer files of the 'tiff' library.")
    sys.exit(1)

setup(
    name='nicoslivewidget',
    version=get_git_version().lstrip('v'),
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
