from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.backends.qt_compat import QtCore, QtWidgets
#import scipy.special._ufuncs_cxx


#_fromUtf8 = QtCore.QString.fromUtf8


try:
    import matplotlib.backends.qt_editor.figureoptions as figureoptions
except ImportError:
    figureoptions = None


class NavigationToolbarPCA(NavigationToolbar2QT, QtWidgets.QToolBar):
    toolitems = (
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to  previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        (None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
        )
    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates)
        self.setIconSize(QtCore.QSize(16, 16))

        self.ct = None
        self.mw = None
        self._idPress1 = None
        self._idPress2 = None
        self._idPress3 = None


    def setController(self, ct):
        self.ct = ct
    def setMainWindow(self, mw):
        self.mw = mw

    def _update_buttons_checked(self):
        # sync button checkstates to match active mode
        self._actions['pan'].setChecked(self._active == 'PAN')
        self._actions['zoom'].setChecked(self._active == 'ZOOM')
        print "updating somehow"
        if self._active is None:

            if self.mw.isCompletePCAView() == False:
                self.exploratoryPCAMode()


    def exploratoryPCAMode(self):
        self.desactiveToolbar()
        print "llamando al pick mode"
        self._idPress = self.canvas.mpl_connect('pick_event', self.ct.onpickpca)
        self._idPress1 = self.canvas.mpl_connect('button_press_event',  self.buttonpress)
        self._idPress2 = self.canvas.mpl_connect('button_release_event',  self.buttonrelease)
        self.mode = 'PCA expl'
        print self.mode
        self.set_message(self.mode)


    def buttonpress(self, event):
        if self._active is None:
            if self._idPress3 is None:
                self._idPress3 = self.canvas.mpl_connect('motion_notify_event',  self.ct.onpickexpl2)

    def buttonrelease(self, event):
        self._idPress3 = self.canvas.mpl_disconnect(self._idPress3)
        self._idPress3 = None


    def desactiveToolbar(self):
        if self.canvas.widgetlock.locked() == True:
            self.canvas.widgetlock.release(self)
        self._idPress = self.canvas.mpl_disconnect(self._idPress)
        self._idPress1 = self.canvas.mpl_disconnect(self._idPress1)
        self._idPress2 = self.canvas.mpl_disconnect(self._idPress2)
        self._active = None
        self._actions['pan'].setChecked(False)
        self._actions['zoom'].setChecked(False)
