# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'FlagConfiguration.ui'
#
# Created: Tue Jul 28 17:05:02 2015
#      by: PyQt4 UI code generator 4.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(400, 300)
        self.groupBox = QtGui.QGroupBox(Dialog)
        self.groupBox.setGeometry(QtCore.QRect(10, 10, 381, 281))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.checkBox = QtGui.QCheckBox(self.groupBox)
        self.checkBox.setGeometry(QtCore.QRect(20, 70, 91, 17))
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.checkBox_2 = QtGui.QCheckBox(self.groupBox)
        self.checkBox_2.setGeometry(QtCore.QRect(20, 40, 91, 17))
        self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        self.lineEdit_6 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_6.setGeometry(QtCore.QRect(300, 40, 51, 20))
        self.lineEdit_6.setObjectName(_fromUtf8("lineEdit_6"))
        self.lineEdit_5 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_5.setGeometry(QtCore.QRect(300, 70, 51, 20))
        self.lineEdit_5.setObjectName(_fromUtf8("lineEdit_5"))
        self.label_8 = QtGui.QLabel(self.groupBox)
        self.label_8.setGeometry(QtCore.QRect(270, 40, 31, 20))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.label_7 = QtGui.QLabel(self.groupBox)
        self.label_7.setGeometry(QtCore.QRect(270, 70, 31, 20))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.lineEdit_3 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_3.setGeometry(QtCore.QRect(180, 40, 51, 20))
        self.lineEdit_3.setObjectName(_fromUtf8("lineEdit_3"))
        self.lineEdit_4 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_4.setGeometry(QtCore.QRect(180, 70, 51, 20))
        self.lineEdit_4.setObjectName(_fromUtf8("lineEdit_4"))
        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setGeometry(QtCore.QRect(140, 40, 51, 20))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setGeometry(QtCore.QRect(140, 70, 51, 20))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.buttonBox = QtGui.QDialogButtonBox(self.groupBox)
        self.buttonBox.setGeometry(QtCore.QRect(20, 240, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.checkBox_3 = QtGui.QCheckBox(self.groupBox)
        self.checkBox_3.setGeometry(QtCore.QRect(20, 100, 70, 17))
        self.checkBox_3.setObjectName(_fromUtf8("checkBox_3"))
        self.checkBox_4 = QtGui.QCheckBox(self.groupBox)
        self.checkBox_4.setGeometry(QtCore.QRect(20, 130, 70, 17))
        self.checkBox_4.setObjectName(_fromUtf8("checkBox_4"))
        self.pushButton = QtGui.QPushButton(self.groupBox)
        self.pushButton.setGeometry(QtCore.QRect(10, 190, 75, 23))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.checkBox_5 = QtGui.QCheckBox(self.groupBox)
        self.checkBox_5.setGeometry(QtCore.QRect(20, 160, 70, 17))
        self.checkBox_5.setObjectName(_fromUtf8("checkBox_5"))
        self.textBrowser = QtGui.QTextBrowser(self.groupBox)
        self.textBrowser.setGeometry(QtCore.QRect(10, 220, 181, 51))
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.groupBox.setTitle(_translate("Dialog", "Threshold settings", None))
        self.checkBox.setText(_translate("Dialog", "PCA threshold", None))
        self.checkBox_2.setText(_translate("Dialog", "Spike rate", None))
        self.lineEdit_6.setText(_translate("Dialog", "100000", None))
        self.lineEdit_5.setText(_translate("Dialog", "100000", None))
        self.label_8.setText(_translate("Dialog", "Max:", None))
        self.label_7.setText(_translate("Dialog", "Max:", None))
        self.lineEdit_3.setText(_translate("Dialog", "0", None))
        self.lineEdit_4.setText(_translate("Dialog", "0", None))
        self.label_6.setText(_translate("Dialog", "Min:", None))
        self.label_5.setText(_translate("Dialog", "Min:", None))
        self.checkBox_3.setText(_translate("Dialog", "Flagged", None))
        self.checkBox_4.setText(_translate("Dialog", "Visible", None))
        self.pushButton.setText(_translate("Dialog", "Clear", None))
        self.checkBox_5.setText(_translate("Dialog", "Selected", None))

