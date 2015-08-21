from PyQt4 import QtGui
from PyQt4.QtGui import QDialog
from matplotlib.backends.qt_compat import QtCore
import numpy as np
import time


class Dialog(QDialog):
    def __init__(self, parent):
        QtGui.QDialog.__init__(self, parent)
        self.setModal(True)

        self.mw = None
        self.ct = None
        self.cw = None

    def clearFields(self):
        self.cw.checkBox.setCheckState(QtCore.Qt.Unchecked)
        self.cw.checkBox_2.setCheckState(QtCore.Qt.Unchecked)
        self.cw.checkBox_3.setCheckState(QtCore.Qt.Unchecked)
        self.cw.checkBox_4.setCheckState(QtCore.Qt.Unchecked)
        self.cw.checkBox_5.setCheckState(QtCore.Qt.Unchecked)



    def setupConfWindow(self, cw, mw, ct):
        self.cw = cw
        self.mw = mw
        self.ct = ct

        QtCore.QObject.connect(self.cw.pushButton, QtCore.SIGNAL('clicked()'), self.clearFields)




    def accept(self):

        inds = np.array(range(len(self.ct.getClusterList())),dtype=int)

        index = self.filterIndexes(inds)
        if index is not None:
            self.ct.filterClusters(index)
            super(Dialog, self).accept()


    def reject(self):
        super(Dialog, self).reject()


    def filterIndexes(self, indexes):

        inds = indexes

        if self.cw.checkBox_2.checkState() == QtCore.Qt.Checked:


            indst = np.where((self.ct.countspikes[inds] > float(self.cw.lineEdit_3.text())) &
                             (self.ct.countspikes[inds] < float(self.cw.lineEdit_6.text())))[0]

            inds = inds[indst]

            if isinstance(inds, int):
                inds = np.array([inds],dtype=int)

        if self.cw.checkBox.checkState() == QtCore.Qt.Checked:

            if self.ct.getPCAList() is None:
                self.cw.textBrowser.setText("PCA not calculated yet")
                return

            indst = np.array([],dtype=int)

            mint = int(self.cw.lineEdit_4.text())
            maxt = int(self.cw.lineEdit_5.text())

            clf = self.ct.getClusterList()


            for i in inds:
                if self.ct.isPCALimits(i, mint, maxt):
                    indst = np.append(indst, i)

            inds = indst


        if self.cw.checkBox_3.checkState() == QtCore.Qt.Checked:

            indst = np.array([],dtype=int)
            for i in self.ct.flagged:
                if i in inds:
                    indst = np.append(indst, i)

            inds = indst


        if self.cw.checkBox_4.checkState() == QtCore.Qt.Checked:


            indp = self.ct.getIndexesVisible(inds)
            inds = inds[indp]

            if isinstance(inds, int):
                inds = np.array([inds],dtype=int)


        if self.cw.checkBox_5.checkState() == QtCore.Qt.Checked:


            indst = np.array([],dtype=int)
            for i in self.ct.getList1Checked():
                if i in inds:
                    indst = np.append(indst, i)

            inds = indst

        return  inds



