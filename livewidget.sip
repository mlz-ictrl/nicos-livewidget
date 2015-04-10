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
//   Georg Brandl <georg.brandl@frm2.tum.de>
//
// *****************************************************************************

%Module(name=nicoslivewidget, version=1)

%Import QtGui/QtGuimod.sip
%Import QtCore/QtCoremod.sip
%Import Qwt5/QwtModule.sip


class LWData
{
%TypeHeaderCode
#include "lw_data.h"
%End
  public:
    LWData(const LWData &other);
    LWData(int width, int height, int depth, const char *data);
    LWData(int width, int height, int depth,
           const char *format, const char *data);
    LWData(const char *filename);

    int width() const;
    int height() const;
    int depth() const;
    double min() const;
    double max() const;

    int currentZ() const;
    virtual void setCurrentZ(int val);

    bool isLog10() const;
    virtual void setLog10(bool val);

    bool hasCustomRange() const;
    double customRangeMin() const;
    double customRangeMax() const;
    virtual void setCustomRange(double lower, double upper);

    virtual double value(double x, double y) const;
    virtual double valueRaw(int x, int y) const;
    virtual double valueRaw(int x, int y, int z) const;

    virtual void histogram(int bins, QVector<double> **xs,
                           QVector<double> **ys) const;
};



class LWZoomer : QwtPlotZoomer
{
%TypeHeaderCode
#include "lw_widget.h"
%End
  public:
    LWZoomer(QwtPlotCanvas *canvas, const QwtPlotSpectrogram *spectro);
};


class LWPlot : QwtPlot
{
%TypeHeaderCode
#include "lw_plot.h"
%End
  public:
    LWPlot(QWidget *parent);

    void updateRange();
    LWZoomer *getZoomer();
    QwtPlotPanner *getPanner();
    QwtPlotPicker *getPicker();
    const QwtRasterData *getData() const;

    void setData(QwtRasterData *data);
    void setGrid(bool val);
    void setColorMap(QwtColorMap &map);

    bool hasGrid();

  public slots:
    void printPlot();
};


class LWWidget : QWidget
{
%TypeHeaderCode
#include "lw_widget.h"
%End
  public:
    LWWidget(QWidget *parent = NULL);
    virtual ~LWWidget();

    void setInstrumentOption(const char *instr);
    int instrument();

    LWPlot *plot();

    LWData *data();
    void setData(LWData *data /Transfer/);

    bool hasGrid() const;
    bool isLog10() const;
    bool isKeepAspect() const;
    bool controlsVisible() const;

    void setCustomRange(double lower, double upper);
    void setStandardColorMap(bool grayscale, bool cyclic);
    void setAxisLabels(const char *xaxis, const char *yaxis);

  protected:
    virtual void resizeEvent(QResizeEvent *event);

  public slots:
    void setGrid(bool val);
    void setLog10(bool val);
    void setKeepAspect(bool val);
    void setControlsVisible(bool val);
    void setControls(LWCtrl which);

    void updateGraph();
    void updateLabels();

  signals:
    void dataUpdated(LWData *data);
    void profilePointPicked(int type, double x, double y);
    void profileUpdate(int type, int nbins, void *dx, void *dy);
};

%UnitCode
#include "lw_common.h"
%End

enum LWFiletype {
    TYPE_NONE,
    TYPE_CASCADE,
    TYPE_FITS,
    TYPE_TOFTOF,
    TYPE_TIFF,
    TYPE_RAW
};

enum LWCtrl {
    Logscale,
    Grayscale,
    Cyclic,
    CreateProfile,
    Histogram,
    MinimumMaximum,
    BrightnessContrast,
    Integrate,
    ImageOperations,
    Despeckle,
    Normalize,
    Darkfield,
    ShowGrid,
    Filelist
};

const char *__version__;

%ModuleCode
#include "version.h"
static const char *__version__ = VERSION;
%End
