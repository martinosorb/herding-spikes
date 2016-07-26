from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.backends.qt_compat import QtCore, QtGui, QtWidgets
import matplotlib
import os
import six
import numpy as np
#import scipy.special._ufuncs_cxx
#_fromUtf8 = QtCore.QString.fromUtf8


try:
    import matplotlib.backends.qt_editor.figureoptions as figureoptions
except ImportError:
    figureoptions = None



class NavigationToolbar3(NavigationToolbar2QT, QtWidgets.QToolBar):
    toolitems = (
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to  previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        (None, None, None, None),
        ('Clusters', 'Select clusters', 'crosshair', 'cursorcluster'),
        ('Selected', 'Shows only the hold on clusters', 'selected', 'selectedmode'),
        (None, None, None, None),
        ('Radio', 'Radial heat map', 'radio', 'radiomode'),
        ('Shape', 'Cluster shape outline', 'shape', 'shapemode'),
        (None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
        )
    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates)
        self.setIconSize(QtCore.QSize(16, 16))

        self.ct = None
        self.mw = None
        self._activeradio = None
        self._activeshape = None
        self._activesel = None
        self._idPress1 = None
        self._idPress2 = None
        self._idPress3 = None




    def setController(self, ct):
        self.ct = ct
    def setMainWindow(self, mw):
        self.mw = mw

    def _init_toolbar(self):
        #self.basedir = os.path.join(matplotlib.rcParams['datapath'], 'images')
        iconDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'icons', '')
        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self.addSeparator()
            else:
                # a = self.addAction(self._icon(iconDir + image_file + '.png'),
                #                          text, getattr(self, callback))
                a = self.addAction(QtGui.QIcon(iconDir + image_file + '.png'),
                                         text, getattr(self, callback))
                self._actions[callback] = a
                if callback in ['zoom', 'pan', 'cursorcluster', 'radiomode', 'shapemode', 'selectedmode']:
                    a.setCheckable(True)
                if tooltip_text is not None:
                    a.setToolTip(tooltip_text)

        if figureoptions is not None:
            # a = self.addAction(self._icon("qt4_editor_options.png"),
            #                    'Customize', self.edit_parameters)
            a = self.addAction(QtGui.QIcon("qt4_editor_options.png"),
                               'Customize', self.edit_parameters)
            a.setToolTip('Edit curves line and axes parameters')

        self.buttons = {}

        # Add the x,y location widget at the right side of the toolbar
        # The stretch factor is 1 which means any resizing of the toolbar
        # will resize this label instead of the buttons.
        if self.coordinates:
            self.locLabel = QtWidgets.QLabel("", self)
            self.locLabel.setAlignment(
                    QtCore.Qt.AlignRight | QtCore.Qt.AlignTop)
            self.locLabel.setSizePolicy(
                QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                  QtWidgets.QSizePolicy.Ignored))
            labelAction = self.addWidget(self.locLabel)
            labelAction.setVisible(True)

        # reference holder for subplots_adjust window
        self.adj_window = None

    if figureoptions is not None:
        def edit_parameters(self):
            allaxes = self.canvas.figure.get_axes()
            if len(allaxes) == 1:
                axes = allaxes[0]
            else:
                titles = []
                for axes in allaxes:
                    title = axes.get_title()
                    ylabel = axes.get_ylabel()
                    label = axes.get_label()
                    if title:
                        fmt = "%(title)s"
                        if ylabel:
                            fmt += ": %(ylabel)s"
                        fmt += " (%(axes_repr)s)"
                    elif ylabel:
                        fmt = "%(axes_repr)s (%(ylabel)s)"
                    elif label:
                        fmt = "%(axes_repr)s (%(label)s)"
                    else:
                        fmt = "%(axes_repr)s"
                    titles.append(fmt % dict(title=title,
                                         ylabel=ylabel, label=label,
                                         axes_repr=repr(axes)))
                item, ok = QtWidgets.QInputDialog.getItem(
                    self.parent, 'Customize', 'Select axes:', titles, 0, False)
                if ok:
                    axes = allaxes[titles.index(six.text_type(item))]
                else:
                    return

            figureoptions.figure_edit(axes, self)


    def cursorcluster(self, *args):
        """Activate the pan/zoom tool. pan with left button, zoom with right"""
        # set the pointer icon and button press funcs to the
        # appropriate callbacks

        if self._active == 'CLUSTER':
            self._active = None
        else:
            self._active = 'CLUSTER'

        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''

        if self._active == 'CLUSTER':
            if self.canvas.widgetlock.locked() == True:
                self.canvas.widgetlock.release(self)
            self._idPress = self.canvas.mpl_connect('pick_event', self.ct.onpickchange)

            self.mode = 'Cluster mode'
            self.canvas.draw()

        self.set_message(self.mode)
        self._update_buttons_checked()


    def radiomode(self, *args):
        """Activate the pan/zoom tool. pan with left button, zoom with right"""
        # set the pointer icon and button press funcs to the
        # appropriate callbacks

        if self._activeradio == 'RADIO':
            self._activeradio = None
        else:
            self._activeradio = 'RADIO'

        if self._activeradio == 'RADIO':
            waslocked = self.canvas.widgetlock.locked()
            if waslocked == True:
                self.canvas.widgetlock.release(self)

            #self.ct.setRadialClustersVisible(True)
            self.mw.heatview(True)
            self._actions['radiomode'].setChecked(True)
            self._actions['selectedmode'].setChecked(False)

            if waslocked == True:
                self.canvas.widgetlock(self)

            self.canvas.draw()
        else:
            waslocked = self.canvas.widgetlock.locked()
            if waslocked == True:
                self.canvas.widgetlock.release(self)



            self.ct.setRadialClustersVisible(False)
            self.mw.clusterview(True)
            self._actions['radiomode'].setChecked(False)

            if waslocked == True:
                self.canvas.widgetlock(self)

            self.canvas.draw()


    def shapemode(self, *args):
        """Activate the pan/zoom tool. pan with left button, zoom with right"""
        # set the pointer icon and button press funcs to the
        # appropriate callbacks

        if self._activeshape == 'SHAPE':
            self._activeshape = None
        else:
            self._activeshape = 'SHAPE'

        if self._activeshape == 'SHAPE':
            waslocked = self.canvas.widgetlock.locked()
            if waslocked == True:
                self.canvas.widgetlock.release(self)

            self.ct.setClusterShapesVisible(True)
            self.ct.enableShapesPicker(True)
            self.ct.enableClustersPicker(False)
            print "checheando shape mode"
            self._actions['shapemode'].setChecked(True)

            if waslocked == True:
                self.canvas.widgetlock(self)

            self.canvas.draw()
        else:
            waslocked = self.canvas.widgetlock.locked()
            if waslocked == True:
                self.canvas.widgetlock.release(self)


            print "unchecheando shape mode"
            self.ct.setClusterShapesVisible(False)
            self.ct.enableShapesPicker(False)
            self.ct.enableClustersPicker(True)
            self._actions['shapemode'].setChecked(False)

            if waslocked == True:
                self.canvas.widgetlock(self)

            self.canvas.draw()


    def selectedmode(self, *args):
        """Activate the pan/zoom tool. pan with left button, zoom with right"""
        # set the pointer icon and button press funcs to the
        # appropriate callbacks

        if self._activesel == 'SELECTED':
            self._activesel = None
        else:
            self._activesel = 'SELECTED'

        if self._activesel == 'SELECTED':
            waslocked = self.canvas.widgetlock.locked()
            if waslocked == True:
                self.canvas.widgetlock.release(self)

            self.ct.initSelectedClusters()

            self._actions['selectedmode'].setChecked(True)


            if waslocked == True:
                self.canvas.widgetlock(self)

            #self.canvas.draw()
        else:
            waslocked = self.canvas.widgetlock.locked()
            if waslocked == True:
                self.canvas.widgetlock.release(self)

            self.ct.finishSelectedClusters()


            self._actions['selectedmode'].setChecked(False)

            if waslocked == True:
                self.canvas.widgetlock(self)

            self.canvas.draw()









    def _update_buttons_checked(self):
        # sync button checkstates to match active mode
        self._actions['pan'].setChecked(self._active == 'PAN')
        self._actions['zoom'].setChecked(self._active == 'ZOOM')
        self._actions['cursorcluster'].setChecked(self._active == 'CLUSTER')

        if self._active is None:
            if self.mw.isClusterView():
                self.exploratoryMode()
            else:
                self.exploratorySpikeMode()

    def desactiveToolbar(self):
        if self.canvas.widgetlock.locked() == True:
            self.canvas.widgetlock.release(self)
        self._idPress = self.canvas.mpl_disconnect(self._idPress)
        self._idPress1 = self.canvas.mpl_disconnect(self._idPress1)
        self._idPress2 = self.canvas.mpl_disconnect(self._idPress2)
        self._active = None
        self._actions['pan'].setChecked(False)
        self._actions['zoom'].setChecked(False)
        self._actions['cursorcluster'].setChecked(False)


    def exploratoryMode(self):
        self.desactiveToolbar()
        self._idPress = self.canvas.mpl_connect('pick_event', self.ct.onpickexpl)
        self.mode = 'Exploratory cluster mode'
        self.set_message(self.mode)




    def exploratorySpikeMode(self):
        self.desactiveToolbar()
        self._idPress = self.canvas.mpl_connect('pick_event', self.ct.onpickspexpl)
        self._idPress1 = self.canvas.mpl_connect('button_press_event',  self.buttonpress)
        self._idPress2 = self.canvas.mpl_connect('button_release_event',  self.buttonrelease)
        self.mode = 'Exploratory spike mode'
        self.set_message(self.mode)

    def buttonpress(self, event):
        if self._active is None:
            if self._idPress3 is None:
                self._idPress3 = self.canvas.mpl_connect('motion_notify_event',  self.ct.onpicksexpl2)

    def buttonrelease(self, event):
        self._idPress3 = self.canvas.mpl_disconnect(self._idPress3)
        self._idPress3 = None


    def isShapeButtonChecked(self):
        return self._actions['shapemode'].isChecked()

    def isSelectButtonChecked(self):
        return self._actions['selectedmode'].isChecked()

    def isHeatButtonChecked(self):
        return self._actions['radiomode'].isChecked()

    def setEnableCursorCluster(self, mode):
        self._actions['cursorcluster'].setEnabled(mode)

    def setEnableClusterShape(self, mode):
        self._actions['shapemode'].setEnabled(mode)

    def setEnableSelectedCluster(self, mode):
        self._actions['selectedmode'].setEnabled(mode)

    def setEnableHeatCluster(self, mode):
        self._actions['radiomode'].setEnabled(mode)

    def isToolbarActive(self):
        active = False
        if self._active is not None:
            active = True
        return active


    def zoom_point(self, xmin, xmax, ymin, ymax):

        # push the current view to define home if stack is empty
        if self._views.empty():
            self.push_current()

        a=self.canvas.figure.get_axes()

        x0 = xmin - 1
        x1 = xmax + 1
        y0 = ymin - 1
        y1 = ymax + 1

        a[0].set_xlim(x0, x1)
        a[0].set_ylim(y0, y1)

        self.draw()
        self.push_current()
