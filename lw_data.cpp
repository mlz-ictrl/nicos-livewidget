// *****************************************************************************
// NICOS, the Networked Instrument Control System of the FRM-II
// Copyright (c) 2009-2014 by the NICOS contributors (see AUTHORS)
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// Module authors:
//   Tobias Weber <tobias.weber@frm2.tum.de>
//   Georg Brandl <georg.brandl@frm2.tum.de>
//   Philipp Schmakat <philipp.schmakat@frm2.tum.de>
//
// *****************************************************************************

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <math.h>
#include <time.h>

#include <fitsio.h>
#include <tiffio.h>

#include <QStringList>

#include "lw_data.h"
#include "lw_imageproc.h"
#include "lw_common.h"

#ifdef CLOCKING
static clock_t clock_start, clock_stop;
#define CLOCK_START()      clock_start = clock()
#define CLOCK_STOP(action) clock_stop  = clock(); \
                           std::cout << (action) << ": " << \
                               (1000 * (float)(clock_stop-clock_start)/CLOCKS_PER_SEC) \
                               << " ms" << std::endl
#else
#define CLOCK_START()
#define CLOCK_STOP(action)
#endif


static inline double safe_log10(double v)
{
    v = (v > 0.) ? log10(v) : -1;
    if (v != v) v = -1.;
    return v;
}

static inline uint16_t bswap_16(uint16_t x)
{
    return (x << 8) | (x >> 8);
}

static inline uint32_t bswap_32(uint32_t x) {
    return (bswap_16(x & 0xFFFF) << 16) | bswap_16(x >> 16);
}

static inline float bswap_32_float(float x)
{
    float res;
    char *d = (char *)&res, *s = (char *)&x;
    *d++ = s[3]; *d++ = s[2]; *d++ = s[1]; *d++ = s[0];
    return res;
}

static inline float bswap_64_float(double x)
{
    double res;
    char *d = (char *)&res, *s = (char *)&x;
    *d++ = s[7]; *d++ = s[6]; *d++ = s[5]; *d++ = s[4];
    *d++ = s[3]; *d++ = s[2]; *d++ = s[1]; *d++ = s[0];
    return res;
}


void replaceExt(std::string& s, const std::string& newExt) {

   std::string::size_type i = s.rfind('.', s.length());

   if (i != std::string::npos) {
      s.replace(i+1, newExt.length(), newExt);
   }
}

/** LWData ********************************************************************/

void LWData::_dummyInit()
{
    m_width = m_height = m_depth = 1;
    // simple constructor: create a 1-element array
    m_data = new data_t[1];
    m_clone = new data_t[1];
    m_data[0] = m_clone[0] = 0;
    m_data_owned = true;
    updateRange();
}

LWData::LWData()
    : m_data(NULL),
      m_clone(NULL),
      m_data_owned(false),
      m_width(1),
      m_height(1),
      m_depth(1),
      m_log10(0),
      m_custom_range(0)
{
    _dummyInit();
}


#define COPY_LOOP(type, data_ptr)                       \
    type *p = (type *)data;                                     \
    data_t *q = data_ptr;                                       \
    for (int i = 0; i < size(); i++){  \
        q[i] = p[i];       \
    }

#define COPY_LOOP_CONVERTED(type, converter, data_ptr)          \
    type *p = (type *)data;                                     \
    data_t *q = data_ptr;                                       \
    for (int i = 0; i < size(); i++) {        \
        q[i] = converter(p[i]);                                 \
    }

void LWData::initFromBuffer(const void *data, std::string format = "<u4")
{
    if (m_data != NULL && m_data_owned) { delete[] m_data;}
    if (m_clone != NULL && m_data_owned) { delete[] m_clone;}
    m_data = new data_t[size()]();
    m_clone = new data_t[size()]();

    m_data_owned = true;
    if (m_data == NULL || m_clone == NULL) {
        std::cerr << "could not allocate memory for data" << std::endl;
        if (m_data != NULL) { delete[] m_data; m_data=NULL;}
        if (m_clone != NULL) { delete[] m_clone;m_clone=NULL;}
        return;
    }
    if (data != NULL) {
    // XXX currently, we interpret signed types as unsigned

    // the easy case
      if (format == "<u4" || format == "<i4" ||
        format == "u4"  || format == "i4") {
        memcpy(m_data, data, sizeof(data_t) * size());
      } else if (format == ">I4" || format == ">i4") {
        COPY_LOOP_CONVERTED(uint32_t, bswap_32, m_data);
      } else if (format == "<u2"  || format == "<i2"  ||
               format == "u2"  || format == "i2" ) {
        COPY_LOOP(uint16_t, m_data);
      } else if (format == "<u1"  || format == "<i1"  ||
               format == "u1"  || format == "i1" ) {
        COPY_LOOP(uint8_t, m_data);
      } else if (format == ">u2" || format == ">i2" ) {
        COPY_LOOP_CONVERTED(uint16_t, bswap_16, m_data);
      } else if (format == "u1" || format == "i1" ) {
        COPY_LOOP(uint8_t, m_data);
      } else if (format == "<f8" || format == "f8" ) {
        COPY_LOOP(double, m_data);
      } else if (format == ">f8" ) {
        COPY_LOOP_CONVERTED(double, bswap_64_float, m_data);
      } else if (format == "<f4" || format == "f4" ) {
        COPY_LOOP(float, m_data);
      } else if (format == ">f4" ) {
        COPY_LOOP_CONVERTED(float, bswap_32_float, m_data);
      } else {
        std::cerr << "Unsupported format: " << format << "!" << std::endl;
      }
    }
    memcpy(m_clone, m_data, sizeof(data_t) * size());
    updateRange();


}


LWData::LWData(int width, int height, int depth, const char *data)
    : m_data(NULL),
      m_clone(NULL),
      m_data_owned(false),
      m_width(width),
      m_height(height),
      m_depth(depth),
      m_cur_z(0),
      m_log10(0),
      m_custom_range(0)
{
    initFromBuffer(data);
}



LWData::LWData(int width, int height, int depth,
               const char *format, const char *data)
    : m_data(NULL),
      m_clone(NULL),
      m_data_owned(false),
      m_width(width),
      m_height(height),
      m_depth(depth),
      m_cur_z(0),
      m_log10(0),
      m_custom_range(0)
{
      initFromBuffer(data, format);

}


LWData::LWData(const char* filename)
    : m_data(NULL),
      m_clone(NULL),
      m_data_owned(false),
      m_width(0),
      m_height(0),
      m_depth(0),
      m_cur_z(0),
      m_log10(0),
      m_custom_range(0),
      m_normalized(0),
      m_darkfieldsubtracted(0),
      m_despeckled(0),
      m_filter(NoImageFilter),
      m_operation(NoImageOperation),
      m_despecklevalue(100),
      m_darkfieldfile(""),
      m_normalizefile("")
{
    if (! _readFits(filename)) {
        if (! _readRaw(filename)) {
           if (! _readTiff(filename)) {
              _dummyInit();
           }
        }
    }
}

LWData::LWData(const LWData &other)
    : m_data(NULL),
      m_clone(NULL),
      m_data_owned(false),
      m_width(other.m_width),
      m_height(other.m_height),
      m_depth(other.m_depth),
      m_min(other.m_min),
      m_max(other.m_max),
      m_cur_z(other.m_cur_z),
      m_log10(other.m_log10),
      m_custom_range(other.m_custom_range),
      m_range_min(other.m_range_min),
      m_range_max(other.m_range_max),
      m_normalized(other.m_normalized),
      m_darkfieldsubtracted(other.m_darkfieldsubtracted),
      m_despeckled(other.m_despeckled),
      m_filter(other.m_filter),
      m_operation(other.m_operation),
      m_despecklevalue(other.m_despecklevalue)
{
    m_data = new data_t[other.size()];
    m_clone = new data_t[other.size()];
    m_data_owned = true;
    memcpy(m_data, other.m_data, sizeof(data_t) * other.size());
    memcpy(m_clone, other.m_clone, sizeof(data_t) * other.size());
}

LWData::~LWData()
{
    if (m_data_owned) {
        if (m_data) {
            delete[] m_data;
            m_data = NULL;
        }
        if (m_clone) {
            delete[] m_clone;
            m_clone = NULL;
        }
        m_data_owned = false;
    }
}

bool LWData::_readFits(const char *filename)
{
    fitsfile *file_pointer;    // CFITSIO file pointer, defined in fitsio.h
    int status = 0;            // CFITSIO status, must be initialized to zero
    int max_dimensions = 3;    // third dimension exists but contains only one image in our case
    int num_dimensions, bitpix, any_null, hdutype;
    float null_value = 0.;

    int total_pixel;
    long dimensions[3];

    float *float_data = NULL;
    data_t *data = NULL;

    if (fits_open_diskfile(&file_pointer, filename, READONLY, &status)) {
        std::cerr << "Could not open file " << filename << " as FITS" <<std::endl;
        return false;
    }
    if (fits_get_img_param(file_pointer, max_dimensions, &bitpix,
                            &num_dimensions, dimensions, &status)) {
        std::cerr << "Could not get image params from " << filename << std::endl;
        fits_close_file(file_pointer, &status);
        return false;
    }
    if (fits_get_hdu_type(file_pointer, &hdutype, &status) ||
        hdutype != IMAGE_HDU || !(num_dimensions == 2 || num_dimensions == 3)) {
        std::cerr << "This .fits file does not contain valid image data!" << std::endl;
        fits_close_file(file_pointer, &status);
        return false;
    }

    m_width  = (int) dimensions[0];
    m_height = (int) dimensions[1];
    m_depth  = 1;

    total_pixel = m_height * m_width;

    float_data = new float[total_pixel]();  // 32bit float values

    if (!float_data) {
        std::cerr << "Memory allocation for data arrays failed!" << std::endl;
        fits_close_file(file_pointer, &status);
        return false;
    }

    if (fits_read_img(file_pointer, TFLOAT, 1, total_pixel, &null_value,
                      float_data, &any_null, &status)) {
        char buf[80];
        fits_read_errmsg(buf);
        std::cerr << "Could not read image data from file: " << buf << std::endl;
        fits_close_file(file_pointer, &status);
        delete[] float_data;
        return false;
    }

    fits_close_file(file_pointer, &status);

    initFromBuffer(float_data, "f4");

    delete[] float_data;
    return true;
}

bool LWData::_readRaw(const char *filename) {

    /* currently assume format is fixed to '<u2'
     * only for the 2-file nicos raw format
     * only 2D-data
    */

    std::ifstream fp;
    std::ifstream hp;
    char* data;
    std::string hname(filename);
    std::string linebuffer;
    size_t totalsize;
    unsigned int dwidth(2);
    std::string type;
    bool success(false);



    replaceExt(hname, "header");
    hp.open(hname.c_str());
    if (!hp) {
       hp.open(filename);
    }
    if (hp) {
        while ( getline(hp,linebuffer)) {
           if (linebuffer.find("ImageType((") != std::string::npos) {
// template for image info
// ImageType((1388, 2064), <type 'numpy.uint16'>, ['X', 'Y'])
// ImageType((1388, 2064), '<u2', ['X', 'Y'])
               char dummy[101];
               sscanf(linebuffer.c_str(), "ImageType((%i, %i), %100s",
                       &m_width, &m_height, dummy);
               int offset=13;
               size_t tstart = linebuffer.find("<type 'numpy.");
               if (tstart == std::string::npos) {
                     tstart = linebuffer.find("'");
                     offset =1;
               }
               size_t tend = linebuffer.find("'", tstart + offset);
               type = linebuffer.substr(tstart + offset, tend - (tstart +offset));
               if (type == "uint16" || type == "int16"  || type == "<u2") {
                   dwidth = 2;
                   type = "<u2";
               }
               if (type == "uint32" || type == "int32" || type == "<u4"){
                  dwidth = 4;
                  type = "<u4";
               }
               m_depth = 1;
               success =true;
               break;
            }
        }
        hp.close();
        if (!success) {
            std::cerr << "Could not read raw header file" << std::endl;
            return false;
        }
    } else {
         std::cerr << "File or header not found" << std::endl;
         return false;
    }
    totalsize = size() * dwidth;
    data = new char[totalsize]();
    if (data!=NULL) {
       fp.open(filename);
       if (!fp) {
           delete[] data;
           std::cerr << "Could not read raw file" << std::endl;
           return false;
       }
       fp.read(data, totalsize);
       if (!fp || fp.gcount() < totalsize) {
           fp.close();
           free(data);
           std::cerr << " Not enough data in raw file" <<std::endl;
           return false;
       }
       fp.close();
       initFromBuffer(data, type);
       delete[] data;
       return true;
    }
    return false;
}

bool LWData::_readTiff(const char* filename) {
    bool success = false;
    TIFF* tif = TIFFOpen(filename, "r");
    if (tif) {
        uint16_t spp, bpp;
        uint32_t linesize;
        char *data;
        std::string itype("<u2");


        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &m_height);
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &m_width);
        if (m_height <1 || m_height > 20000 || m_width <1 || m_width > 20000) {
           // limit to at least 1x1 to max 20000x20000 pixel
           TIFFClose(tif);
           return false;
        }
        m_depth = 1;
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bpp);
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
        if (bpp == 16) {
           itype = "<u2";
        } else if (bpp == 8 ) {
           itype = "<u1";
        }
        linesize = TIFFScanlineSize(tif);
        data = new char[linesize * m_height]();
        if (!data) {
           TIFFClose(tif);
           return false;
        }

        for (size_t i = 0; i < m_height; i++) {
               TIFFReadScanline(tif, &data[i * linesize], i, 0);
        }

        TIFFClose(tif);
        initFromBuffer(data, itype);
        delete[] data;
        success = true;
   }
   return success;
}



inline data_t LWData::data(int x, int y, int z) const
{
    if (m_data == NULL)
        return 0;
    if (x >= 0 && x < m_width &&
        y >= 0 && y < m_height &&
        z >= 0 && z < m_depth)
        return m_data[z*m_width*m_height + y*m_width + x];
    return 0;
}

void LWData::updateRange()
{
    m_min = std::numeric_limits<double>::max();
    m_max = 0;
    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            double v = value((double)x, (double)y);
            m_min = (m_min < v) ? m_min : v;
            m_max = (m_max > v) ? m_max : v;
        }
    }
}

double LWData::value(double x, double y) const
{
    double v = (double)data((int)x, (int)y, m_cur_z);
    if (m_log10)
        v = safe_log10(v);
    /*
    if (m_custom_range) {
        v = (v > m_range_max) ? m_range_max : v;
        v = (v < m_range_min) ? m_range_min : v;
    }
    */
    return v;
}

double LWData::valueRaw(int x, int y) const
{
    return (double)data(x, y, m_cur_z);
}

double LWData::valueRaw(int x, int y, int z) const
{
    return (double)data(x, y, z);
}

void LWData::histogram(int bins, double *xs, double *ys) const
{
    double step = (m_max - m_min) / (double)bins;
    if (step == 0)
        return;
    for (int i = 0; i < bins; ++i) {
        xs[i] = m_min + i * step + 0.5 * step;
    }
    std::fill(ys, ys+bins, 0.0);
    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            ys[(int)((value(x, y) - m_min) / step)]++;
        }
    }
}

void LWData::histogram(int bins, QVector<double> **xs, QVector<double> **ys) const
{
    *xs = new QVector<double>(bins);
    *ys = new QVector<double>(bins + 1);
    double step = (m_max - m_min) / (double)bins;
    if (step == 0)
        return;
    for (int i = 0; i < bins; ++i) {
        (**xs)[i] = m_min + i * step + 0.5 * step;
    }
    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            (**ys)[(int)((value(x, y) - m_min) / step)]++;
        }
    }
    (**ys)[bins-1] += (**ys)[bins];
    (*ys)->pop_back();
}

void LWData::setCurrentZ(int val)
{
    if (val < 0 || val >= m_depth) {
        std::cerr << "invalid current Z selected" << std::endl;
        return;
    }
    m_cur_z = val;
    updateRange();
}

void LWData::setLog10(bool val)
{
    if (m_log10 != val) {
        if (m_custom_range) {
            if (val) {
                m_range_min = safe_log10(m_range_min);
                m_range_max = safe_log10(m_range_max);
            } else {
                m_range_min = exp(m_range_min * log(10.));
                m_range_max = exp(m_range_max * log(10.));
            }
        }
        m_log10 = val;
        updateRange();
    }
}

void LWData::setDespeckled(bool val)
{
    if (m_despeckled == val)
        return;
    m_despeckled = val;

    if (m_despeckled) {
        CLOCK_START();
        float *pdata = (float *)malloc(size() * sizeof(float));
        for (int i = 0; i < size(); ++i)
            pdata[i] = (float)m_data[i];
        CLOCK_STOP("malloc and copy");

        CLOCK_START();
        LWImageProc::despeckleFilter(pdata, m_despecklevalue, m_width, m_height);
        CLOCK_STOP("despeckle filter");

        CLOCK_START();
        for (int i = 0; i < size(); ++i)
            m_data[i] = (data_t)pdata[i];
        CLOCK_STOP("copy back to data");

        free(pdata);
    } else {
        memcpy(m_data, m_clone, sizeof(data_t) * size());
    }
    updateRange();
}

void LWData::setDespeckleValue(float value)
{
    if (m_despecklevalue == value)
        return;
    m_despecklevalue = value;

    memcpy(m_data, m_clone, sizeof(data_t) * size());

    float *pdata = (float *)malloc(size() * sizeof(float));
    for (int i = 0; i < size(); ++i)
        pdata[i] = (float)m_data[i];

    LWImageProc::despeckleFilter(pdata, m_despecklevalue, m_width, m_height);

    for (int i = 0; i < size(); ++i)
        m_data[i] = (data_t)pdata[i];

    free(pdata);

    updateRange();
}

void LWData::setNormalized(bool val)
{
    if (m_normalized == val)
        return;
    m_normalized = val;

    if (m_normalized) {
        float *data = (float *)malloc(size() * sizeof(float));
        for (int i = 0; i < size(); ++i)
            data[i] = (float)m_data[i];

        CLOCK_START();
        LWData openbeam(m_normalizefile.toStdString().c_str());
        float *ob_data = (float *)malloc(size() * sizeof(float));
        for (int i = 0; i < size(); ++i)
            ob_data[i] = openbeam.buffer()[i];
        CLOCK_STOP("loaded openbeam image");

        CLOCK_START();
        LWData darkfield(m_darkfieldfile.toStdString().c_str());
        float *di_data = (float *)malloc(size() * sizeof(float));
        for (int i = 0; i < size(); ++i)
            di_data[i] = darkfield.buffer()[i];
        CLOCK_STOP("loaded dark image");

        if (m_despeckled) {
            CLOCK_START();
            LWImageProc::despeckleFilter(di_data, m_despecklevalue, m_width, m_height);
            LWImageProc::despeckleFilter(ob_data, m_despecklevalue, m_width, m_height);
            LWImageProc::despeckleFilter(data, m_despecklevalue, m_width, m_height);
            CLOCK_STOP("removed gamma spots");
        }

        CLOCK_START();
        LWImageProc::pixelwiseSubtractImages(ob_data, di_data, m_width, m_height);
        CLOCK_STOP("pixelwise subtract dark image from openbeam image");

        CLOCK_START();
        LWImageProc::pixelwiseSubtractImages(data, di_data, m_width, m_height);
        CLOCK_STOP("pixelwise subtract dark image from data");

        CLOCK_START();
        LWImageProc::pixelwiseDivideImages(data, ob_data, m_width, m_height);
        CLOCK_STOP("pixelwise divide images");

        clampedCopyFloatVals(data);

        free(data);
        free(ob_data);
        free(di_data);

        updateRange();
    } else {
        CLOCK_START();
        memcpy(m_data, m_clone, sizeof(data_t) * size());
        CLOCK_STOP("restore original from memory");

        CLOCK_START();
        updateRange();
        CLOCK_STOP("updateRange() function");
    }
}


void LWData::setNormalizeFile(QString val)
{
    if (m_normalizefile == val)
        return;
    m_normalizefile = val;
}


void LWData::clampedCopyFloatVals(float* pdata){
        for (int i = 0; i < size(); ++i){
            if ( pdata[i] >  std::numeric_limits<data_t>::min()) {
		if ( pdata[i] < std::numeric_limits<data_t>::max() ) {
	            m_data[i] = (data_t)pdata[i];
                } else {
                    m_data[i] = std::numeric_limits<data_t>::max() ;
                }
            } else {
                m_data[i] = std::numeric_limits<data_t>::min() ;
           }
        }
}

void LWData::setDarkfieldSubtracted(bool val)
{
    if (m_darkfieldsubtracted == val)
        return;
    m_darkfieldsubtracted = val;

    if (m_darkfieldsubtracted) {
        float *pdata = (float *)malloc(size() * sizeof(float));
        for (int i = 0; i < size(); ++i)
            pdata[i] = (float)m_data[i];

        CLOCK_START();
        LWData darkfield(m_darkfieldfile.toStdString().c_str());
        float *sdata = (float *)malloc(size() * sizeof(float));
        for (int i = 0; i < size(); ++i)
            sdata[i] = darkfield.buffer()[i];
        CLOCK_STOP("load darkfield image");

        CLOCK_START();
        LWImageProc::pixelwiseSubtractImages(pdata, sdata, m_width, m_height);
        CLOCK_STOP("pixelwise subtract images");

        clampedCopyFloatVals(pdata);
        free(pdata);
        free(sdata);
        updateRange();
    } else {
        CLOCK_START();
        memcpy(m_data, m_clone, sizeof(data_t) * size());
        CLOCK_STOP("restore original from memory");

        CLOCK_START();
        updateRange();
        CLOCK_STOP("updateRange() function");
    }
}

void LWData::setDarkfieldFile(QString val)
{
    if (m_darkfieldfile == val)
        return;
    m_darkfieldfile = val;
}

void LWData::setImageFilter(LWImageFilters which)
{
    if (m_filter == which)
        return;

    m_filter = which;
    if (m_filter == NoImageFilter) {
        memcpy(m_data, m_clone, sizeof(data_t) * size());
    } else {
        float *pdata = (float *)malloc(size() * sizeof(float));

        for (int i=0; i<size(); ++i)
            pdata[i] = (float)m_data[i];

        if (m_filter == MedianFilter) {
            LWImageProc::medianFilter(pdata, m_width, m_height);
        } else if (m_filter == HybridMedianFilter) {
            LWImageProc::hybridmedianFilter(pdata, m_width, m_height);
        } else if (m_filter == DespeckleFilter) {
            LWImageProc::despeckleFilter(pdata, m_despecklevalue, m_width, m_height);
        }

        for (int i = 0; i < size(); ++i)
            m_data[i] = (data_t)pdata[i];

        free(pdata);
    }
    updateRange();
}

void LWData::setImageOperation(LWImageOperations which)
{
    if (m_operation == which)
        return;

    m_operation = which;
    if (m_operation == NoImageOperation) {
        memcpy(m_data, m_clone, sizeof(data_t) * size());
    } else {
        float *pdata = (float *)malloc(size() * sizeof(float));

        for (int i = 0; i < size(); ++i)
            pdata[i] = (float)m_data[i];

        if (m_operation == StackAverage) {
            str_vec myList;

            LWImageProc::pixelwiseAverage(pdata, myList, m_width, m_height);
        }

        for (int i = 0; i < size(); ++i)
            m_data[i] = (data_t)pdata[i];

        free(pdata);
    }
    updateRange();
}


double LWData::customRangeMin() const
{
    if (m_custom_range)
        return m_range_min;
    return m_min;
}

double LWData::customRangeMax() const
{
    if (m_custom_range)
        return m_range_max;
    return m_max;
}

void LWData::setCustomRange(double lower, double upper)
{
    if (lower == 0 && upper == 0) {
        m_custom_range = false;
    } else {
        m_custom_range = true;
        m_range_min = (lower < upper) ? lower : upper;
        m_range_max = (lower < upper) ? upper : lower;
    }
    updateRange();
}



void LWData::saveAsFitsImage(float *data, char *fits_filename)
{
    long 		 dimensions[2] = {m_width, m_height};
    fitsfile    *ffptr;
    int 		 status = 0;
    int          exists = 0;

    if (fits_file_exists(fits_filename, &exists, &status)) {
        if (fits_create_file(&ffptr, fits_filename, &status))
            std::cerr << "doesn't exist but still couldn't create file" << std::endl;
    } else {
        // XXX doesn't belong here, LWData is not supposed to use GUI
/*        switch (QMessageBox::warning(NULL, "Replace file?",
                                     "Do you want to replace \n" + (QString)fits_filename, "&Yes", "&No") )
        {
            case 0:
                if (fits_open_file(&ffptr, fits_filename, READWRITE, &status))
                    std::cerr << "couldn't open file" << std::endl;
                if (fits_delete_file(ffptr, &status))
                    std::cerr << "couldn't delete existing file" << std::endl;
                if (fits_create_file(&ffptr, fits_filename, &status))
                    std::cerr << "couldn't create file" << std::endl;
                break;

            default:
                break;
        }
*/
    }

    if (fits_create_img(ffptr, 32, 2, dimensions, &status))
        std::cerr << "couldn't add image table" << std::endl;

    if (fits_write_img(ffptr, TFLOAT, 1, 2048*2048, data, &status))
        std::cerr << "couldn't write data" << std::endl;

    fits_close_file(ffptr, &status);

    if (status)
        fits_report_error(stderr, status);   // print any cfitsio error message
}



float LWData::getFloatFromFitsHeader(const char *filename, const char *headerEntry)
{
    fitsfile *file_pointer;    // CFITSIO file pointer, defined in fitsio.h
    int status = 0;            // CFITSIO status, must be initialized to zero

    char value[FLEN_CARD];

    if (!fits_open_file(&file_pointer, filename, READONLY, &status))
    {
        if (!fits_read_card(file_pointer, headerEntry, value, &status))
        {
            fits_close_file(file_pointer, &status);

            if (status)
                fits_report_error(stderr, status);   // print any cfitsio error message

            std::cout << value << std::endl;
            QString val = value;
            QStringList line = val.split(" ");
            float v;

            if (line.at(0).toStdString() != "HIERARCH")
            {
                line = val.split("=");
                v = line.at(1).toFloat();
                std::cout << line.at(1).toFloat() << std::endl;
            }
            else
            {
                line = val.split("'");
                v = line.at(1).split(" ")[0].toFloat();
                std::cout << line.at(1).split(" ")[0].toFloat() << std::endl;
            }

            return v;
        }
        else
        {
            std::cout << "Could not find header keyword: " << headerEntry << std::endl;
            return 0.0;
        }
    }
    else
    {
        std::cerr << "Could not open file " << filename << std::endl;
        return 0.0;
    }
}

std::string LWData::getStringFromFitsHeader(const char *filename, const char *headerEntry)
{
    fitsfile *file_pointer;    // CFITSIO file pointer, defined in fitsio.h
    int status = 0;            // CFITSIO status, must be initialized to zero

    char value[FLEN_CARD];

    if (!fits_open_file(&file_pointer, filename, READONLY, &status))
    {
        if (!fits_read_card(file_pointer, headerEntry, value, &status))
        {
            fits_close_file(file_pointer, &status);

            if (status)
                fits_report_error(stderr, status);   // print any cfitsio error message

            std::cout << value << std::endl;
            QString val = value;
            QStringList line = val.split(" ");
            std::string v;

            if (line.at(0).toStdString() == "HIERARCH")
            {
                line = val.split("'");
                v = line.at(1).split(" ")[0].toStdString();
                std::cout << v << std::endl;
            }
            else
            {
                v = "+";
                std::cout << "default polarity is positive" << std::endl;
            }

            return v;
        }
        else
        {
            std::cout << "Could not find header keyword: " << headerEntry << std::endl;
            return " ";
        }
    }
    else
    {
        std::cerr << "Could not open file " << filename << std::endl;
        return " ";
    }
}
