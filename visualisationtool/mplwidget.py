import string
from PyQt4 import QtGui
from PyQt4.QtCore import Qt, QAbstractTableModel, QVariant, QAbstractItemModel
from PyQt4.QtGui import QColor, QPalette, QLabel, QPixmap
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from PyQt4 import QtCore
from customtoolbar import NavigationToolbar3
from customtoolbarpca import NavigationToolbarPCA
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
#import scipy.special._ufuncs_cxx
from matplotlib.figure import Figure

#_fromUtf8 = QtCore.QString.fromUtf8


class MplCanvas(FigureCanvas):
    def __init__(self):
        # self.fig = plt.figure(figsize=(20, 20))
        # plt.ion()
        self.fig = Figure(figsize=(20,20))#, dpi=100)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                    QtGui.QSizePolicy.Expanding,
                                    QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.ax = self.fig.add_subplot(111)
        # self.ax.hold(False)


class MplWidget(QtGui.QWidget):
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.ntb = NavigationToolbar3(self.canvas, self)
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.ntb)
        self.setLayout(self.vbl)

class MplWidgetT(QtGui.QWidget):
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.ntb = NavigationToolbar(self.canvas, self)
        self.ntb.setIconSize(QtCore.QSize(16, 16))
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.ntb)
        self.setLayout(self.vbl)

class MplWidgetPCA(QtGui.QWidget):
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.ntb = NavigationToolbarPCA(self.canvas, self)
        self.ntb.setIconSize(QtCore.QSize(16, 16))
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.ntb)
        self.setLayout(self.vbl)

class MplWidget2(QtGui.QWidget):
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self,parent)
        self.canvas = MplCanvas()
        self.checkBox = QtGui.QCheckBox()
        self.checkBox.setGeometry(QtCore.QRect(40, 30, 70, 17))
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.checkBox)
        self.setLayout(self.vbl)



class MplWidget3(QtGui.QWidget):
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)


class WidgetItem1(QtGui.QListWidgetItem):
    def __init__(self, parent=None):
        QtGui.QListWidgetItem.__init__(self,parent)
        self.setFlags(self.flags() | QtCore.Qt.ItemIsUserCheckable)
        self.setCheckState(QtCore.Qt.Unchecked)


    def setColour(self, color):
        self.setBackgroundColor(QColor(int(color[0]*255), int(color[1]*255.0), int(color[2]*255.0), 100))



class WidgetCheckBox1(QtGui.QLabel):
    def __init__ (self, parent = None):
        super(WidgetCheckBox1, self).__init__(parent)

        pix= QPixmap(_fromUtf8('checkbox.png'))
        self.state = QtCore.Qt.Unchecked
        self.setPixmap(pix)
        self.setAlignment(Qt.AlignCenter)

    def changeState(self):
        if self.state == QtCore.Qt.Unchecked:
            pix= QPixmap(_fromUtf8('checkboxcheck.png'))
            self.state = QtCore.Qt.Checked
        else:
            pix= QPixmap(_fromUtf8('checkbox.png'))
            self.state = QtCore.Qt.Unchecked

        self.setPixmap(pix)
        self.setAlignment(Qt.AlignCenter)

    def setCheckState(self, state):
        if state == QtCore.Qt.Unchecked:
            pix= QPixmap(_fromUtf8('checkbox.png'))
            self.state = QtCore.Qt.Unchecked
        else:
            pix= QPixmap(_fromUtf8('checkboxcheck.png'))
            self.state = QtCore.Qt.Checked

        self.setPixmap(pix)
        self.setAlignment(Qt.AlignCenter)

    def checkState(self):
        return self.state



class WidgetFlag(QtGui.QLabel):
    def __init__ (self, parent = None):
        super(WidgetFlag, self).__init__(parent)
        pix= QPixmap(_fromUtf8('empty_flag.png'))
        self.setPixmap(pix)
        self.state = QtCore.Qt.Unchecked
        self.setAlignment(Qt.AlignCenter)


    def changeState(self):
        if self.state == QtCore.Qt.Unchecked:
            pix= QPixmap(_fromUtf8('flag_checked.png'))
            self.state = QtCore.Qt.Checked
        else:
            pix= QPixmap(_fromUtf8('empty_flag.png'))
            self.state = QtCore.Qt.Unchecked

        self.setPixmap(pix)
        self.setAlignment(Qt.AlignCenter)

    def setCheckState(self, state):
        if state == QtCore.Qt.Unchecked:
            pix= QPixmap(_fromUtf8('empty_flag.png'))
            self.state = QtCore.Qt.Unchecked
        else:
            pix= QPixmap(_fromUtf8('flag_checked.png'))
            self.state = QtCore.Qt.Checked

        self.setPixmap(pix)
        self.setAlignment(Qt.AlignCenter)

    def checkState(self):
        return self.state











class WidgetItem2(QtGui.QListWidgetItem):
    def __init__(self, parent=None):
        QtGui.QListWidgetItem.__init__(self,parent)
        self.setFlags(self.flags() | QtCore.Qt.ItemIsUserCheckable)
        self.setCheckState(QtCore.Qt.Unchecked)
        #self.setSizeHint(QtCore.QSize(120,60))
        #self.setHidden(True)

    def setColour(self, color):
        self.setBackgroundColor(QColor(int(color[0]*255), int(color[1]*255.0), int(color[2]*255.0), 100))

class TableWidgetItem2(QtGui.QTableWidgetItem):

    def __init__(self, color1 = 'white', color2 = 'black', text1 = '', text2 = '',parent=None):
        QtGui.QTableWidgetItem.__init__(self)
        self.setBackgroundColor(QColor(color1))
        self.state = QtCore.Qt.Unchecked
        self.color1 =color1
        self.color2 =color2
        self.textSel = text2
        self.textUn = text1
        self.setText(self.textUn)


    def changeState(self):
        if self.state == QtCore.Qt.Unchecked:
            self.setBackgroundColor(QColor(self.color2))
            self.setText(self.textSel)
            self.state = QtCore.Qt.Checked
        else:
            self.setBackgroundColor(QColor(self.color1))
            self.setText(self.textUn)
            self.state = QtCore.Qt.Unchecked

    def setCheckState(self, state):
        if state == QtCore.Qt.Unchecked:
            self.setBackgroundColor(QColor(self.color1))
            self.setText(self.textUn)
            self.state = QtCore.Qt.Unchecked
        else:
            self.setBackgroundColor(QColor(self.color2))
            self.setText(self.textSel)
            self.state = QtCore.Qt.Checked

    def checkState(self):
        return self.state

    def __lt__(self, other):
        return False


class TableWidgetItem1(QtGui.QTableWidgetItem):

    def __init__(self, text = '',parent=None):
        QtGui.QTableWidgetItem.__init__(self, text)

    def __lt__(self, other):
        #print "data"
        #print self.text().toLongLong()
        return (float(self.text().replace("*", "")) <
                float(other.text().replace("*", "")))

class TableWidgetItem4(QtGui.QTableWidgetItem):

    def __init__(self, text = '',parent=None):
        QtGui.QTableWidgetItem.__init__(self, text)

    def __lt__(self, other):

        return (float(self.text()) <
                float(other.text()))


class TableWidgetItem3(QtGui.QTableWidgetItem):

    def __init__(self, color1 = 'white', color2 = 'black',parent=None):
        QtGui.QTableWidgetItem.__init__(self)
        self.setBackgroundColor(QColor(color1))
        self.state = QtCore.Qt.Unchecked
        self.color1 =color1
        self.color2 =color2


    def changeState(self):
        if self.state == QtCore.Qt.Unchecked:
            self.setBackgroundColor(QColor(self.color2))
            self.state = QtCore.Qt.Checked
        else:
            self.setBackgroundColor(QColor(self.color1))
            self.state = QtCore.Qt.Unchecked

    def setCheckState(self, state):
        if state == QtCore.Qt.Unchecked:
            self.setBackgroundColor(QColor(self.color1))
            self.state = QtCore.Qt.Unchecked
        else:
            self.setBackgroundColor(QColor(self.color2))
            self.state = QtCore.Qt.Checked

    def checkState(self):
        return self.state

    def __lt__(self, other):
        return False


class MySortFilterProxyModel(QtGui.QSortFilterProxyModel):
    def __init__(self, parent=None):
        super(MySortFilterProxyModel, self).__init__(parent)
        self.rows = []

    def setHHeader(self, header):
        self.hheader = header

    def setSource(self, source):
        self.source = source


    #def lessThan(self, QModelIndex, QModelIndex_1):
    #    return float(self.source.data(QModelIndex, Qt.DisplayRole)) < \
    #           float(self.source.data(QModelIndex_1, Qt.DisplayRole))






    def headerData(self, p_int, Qt_Orientation, int_role=None):
        if(Qt_Orientation == Qt.Horizontal):
            if(int_role == Qt.DisplayRole):
                return self.hheader[p_int]




        if(Qt_Orientation == Qt.Vertical):
            if(int_role == Qt.DisplayRole):
                return p_int

    def setValidRows(self, rows):
        self.rows = rows

    def getValidRows(self):
        return self.rows

    def filterAcceptsRow(self, sourceRow, sourceParent):
        return self.rows[sourceRow]




class MyTableModel(QAbstractTableModel):
    def __init__(self, datain, headerdata, parent=None, *args):
        """ datain: a list of lists
            headerdata: a list of strings
        """
        QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = datain
        self.headerdata = headerdata
        self.holdonColor = 'grey'
        self.holdOffColor= 'white'
        self.flagColor = 'orange'
        self.unFlagColor = 'green'
        self.selectColor = 'blue'
        self.selstatus = ['U', 'H', 'U', 'H']
        self.selcol  = [self.holdOffColor, self.holdonColor, self.selectColor, self.selectColor]
        self.flagstatus = ['U', 'F']
        self.flagcol = [self.unFlagColor, self.flagColor]
        self.sortstring = [0,3,4]
        self.sortnumber = [1,2]
        self.tableLength = 0

    def setPCAValid(self):
        self.sortstring = [0,4]
        self.sortnumber = [1,2,3]

    def flags(self, index):

        flags = QAbstractItemModel.flags(self, index) ^ QtCore.Qt.ItemIsSelectable
        return flags



    def setColors(self, colors):
        self.color = colors
        self.flagcolors = [self.unFlagColor]*len(colors)
        self.holdcolors = [self.holdOffColor]*len(colors)



    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        return len(self.arraydata[0])

    ##def setTableLength(self, length):
    ##    self.tableLength = length

    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        #if (index.row() >= len(self.arraydata))
        elif (role == Qt.DisplayRole) | (role == Qt.EditRole):
            return self.arraydata[index.row()][index.column()]


        elif (role == Qt.BackgroundRole):
            if index.column() == 0:
                return QColor(self.holdcolors[int(str(self.arraydata[index.row()][1]))])

            if index.column() == 4:
                return QColor(self.flagcolors[int(str(self.arraydata[index.row()][1]))])

            if index.column() == 1:
                cl = self.arraydata[index.row()][1]
                color = self.color[int(str(cl))]
                return QColor(int(color[0]*255), int(color[1]*255.0), int(color[2]*255.0), 100)

        elif (role == Qt.UserRole):

            if index.column() in self.sortstring:
                #return self.arraydata[index.row()][index.column()]
                return 0

            elif index.column() in self.sortnumber:
                    return float(self.arraydata[index.row()][index.column()])





        return QVariant()

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if role == QtCore.Qt.EditRole:
            val = value.toString()

            if index.column() == 0:
                self.arraydata[index.row()][0] = self.selstatus[int(val)]
                self.holdcolors[int(self.arraydata[index.row()][1])] = self.selcol[int(val)]

            if index.column() == 3:
                self.arraydata[index.row()][3] = val


            if index.column() == 4:
                self.arraydata[index.row()][4] = self.flagstatus[int(val)]
                self.flagcolors[int(self.arraydata[index.row()][1])] = self.flagcol[int(val)]



            self.dataChanged.emit(index, index)
            return True




        return False

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headerdata[col])
        return QVariant()
