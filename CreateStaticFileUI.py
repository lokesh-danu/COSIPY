# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'CreateStaticFile.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(446, 259)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(26, 19, 401, 171))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayoutWidget = QtWidgets.QWidget(self.tab)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(0, 10, 391, 121))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setVerticalSpacing(10)
        self.gridLayout.setObjectName("gridLayout")
        self.TOutputName = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.TOutputName.setObjectName("TOutputName")
        self.gridLayout.addWidget(self.TOutputName, 3, 1, 1, 1)
        self.BDem = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.BDem.setObjectName("BDem")
        self.gridLayout.addWidget(self.BDem, 0, 2, 1, 1)
        self.TShape = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.TShape.setObjectName("TShape")
        self.gridLayout.addWidget(self.TShape, 1, 1, 1, 1)
        self.TDem = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.TDem.setObjectName("TDem")
        self.gridLayout.addWidget(self.TDem, 0, 1, 1, 1)
        self.BShape = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.BShape.setObjectName("BShape")
        self.gridLayout.addWidget(self.BShape, 1, 2, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 3, 0, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_9.setObjectName("label_9")
        self.gridLayout.addWidget(self.label_9, 2, 0, 1, 1)
        self.BOutput = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.BOutput.setObjectName("BOutput")
        self.gridLayout.addWidget(self.BOutput, 2, 2, 1, 1)
        self.TOutput = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.TOutput.setObjectName("TOutput")
        self.gridLayout.addWidget(self.TOutput, 2, 1, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.tab_2)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(0, 10, 211, 123))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_3.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_4 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_4.setObjectName("label_4")
        self.gridLayout_3.addWidget(self.label_4, 1, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_6.setObjectName("label_6")
        self.gridLayout_3.addWidget(self.label_6, 3, 0, 1, 1)
        self.TUlat = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.TUlat.setEnabled(True)
        self.TUlat.setObjectName("TUlat")
        self.gridLayout_3.addWidget(self.TUlat, 1, 1, 1, 1)
        self.TLLat = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.TLLat.setEnabled(True)
        self.TLLat.setObjectName("TLLat")
        self.gridLayout_3.addWidget(self.TLLat, 2, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_5.setObjectName("label_5")
        self.gridLayout_3.addWidget(self.label_5, 2, 0, 1, 1)
        self.TLLon = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.TLLon.setEnabled(True)
        self.TLLon.setObjectName("TLLon")
        self.gridLayout_3.addWidget(self.TLLon, 3, 1, 1, 1)
        self.TRLon = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.TRLon.setEnabled(True)
        self.TRLon.setObjectName("TRLon")
        self.gridLayout_3.addWidget(self.TRLon, 4, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_7.setObjectName("label_7")
        self.gridLayout_3.addWidget(self.label_7, 4, 0, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_10.setObjectName("label_10")
        self.gridLayout_3.addWidget(self.label_10, 0, 0, 1, 2)
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.tab_2)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(220, 10, 160, 51))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout_4.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_8 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_8.setObjectName("label_8")
        self.gridLayout_4.addWidget(self.label_8, 1, 0, 1, 1)
        self.TAggregateDegree = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.TAggregateDegree.setEnabled(True)
        self.TAggregateDegree.setObjectName("TAggregateDegree")
        self.gridLayout_4.addWidget(self.TAggregateDegree, 1, 1, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_11.setObjectName("label_11")
        self.gridLayout_4.addWidget(self.label_11, 0, 0, 1, 2)
        self.adjustVals = QtWidgets.QPushButton(self.tab_2)
        self.adjustVals.setGeometry(QtCore.QRect(260, 80, 75, 23))
        self.adjustVals.setObjectName("adjustVals")
        self.tabWidget.addTab(self.tab_2, "")
        self.CreateStaticFile = QtWidgets.QPushButton(self.centralwidget)
        self.CreateStaticFile.setGeometry(QtCore.QRect(140, 200, 161, 41))
        self.CreateStaticFile.setObjectName("CreateStaticFile")
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Create Static File"))
        self.BDem.setText(_translate("MainWindow", "Browse"))
        self.BShape.setText(_translate("MainWindow", "Browse"))
        self.label_2.setText(_translate("MainWindow", "DEM Path"))
        self.label.setText(_translate("MainWindow", "Shape File Path"))
        self.label_3.setText(_translate("MainWindow", "Output Static Filename"))
        self.label_9.setText(_translate("MainWindow", "Output Static Folder Path"))
        self.BOutput.setText(_translate("MainWindow", "Browse"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Set Paths"))
        self.label_4.setText(_translate("MainWindow", "Upper Latitude"))
        self.label_6.setText(_translate("MainWindow", "Left Longitude"))
        self.label_5.setText(_translate("MainWindow", "Lower Latitude"))
        self.label_7.setText(_translate("MainWindow", "Right Longitude"))
        self.label_10.setText(_translate("MainWindow", "Shrink DEM ( Lat/Lon Limits of ERA 5 Data)"))
        self.label_8.setText(_translate("MainWindow", "Aggregate Degree"))
        self.label_11.setText(_translate("MainWindow", "            Aggregate DEM"))
        self.adjustVals.setText(_translate("MainWindow", "Adjust Values"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Settings"))
        self.CreateStaticFile.setText(_translate("MainWindow", "Create Static File"))
