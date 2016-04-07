#!/usr/bin/python

import sys
from WindowDesign import *
from MainWindow import MainWindow


if __name__ == "__main__":
	app =QtGui.QApplication(sys.argv)
	myapp = MainWindow()
	myapp.show()
	sys.exit(app.exec_())
