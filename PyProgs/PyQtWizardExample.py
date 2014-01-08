__author__ = 'amarch'


from PyQt4 import QtCore
from PyQt4 import QtGui
from PyQt4.QtCore import pyqtProperty

class QIComboBox(QtGui.QComboBox):
    def __init__(self,parent=None):
        super(QIComboBox, self).__init__(parent)

    @property
    def currentItemData(self):
        return self.itemData(self.currentIndex()).toString()

class VariantWizard(QtGui.QWizard):
    def __init__(self, parent=None):
        super(VariantWizard, self).__init__(parent)
        self.addPage(Page1(self))
        self.addPage(Page2(self))
        self.setWindowTitle("QVariant Test")
        self.resize(640,480)

class Page1(QtGui.QWizardPage):
    def __init__(self, parent=None):
        super(Page1, self).__init__(parent)
        self.version_combo = QIComboBox(self)
        self.version_combo.addItem("filename1","/path/to/filename1")
        self.version_combo.addItem("filename2","/path/to/filename2")
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.version_combo)
        self.setLayout(layout)

        self.registerField("version",self.version_combo, "currentItemData")

class Page2(QtGui.QWizardPage):
    def __init__(self, parent=None):
        super(Page2, self).__init__(parent)
        self.label1 = QtGui.QLabel()
        self.label2 = QtGui.QLabel()
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.label1)
        layout.addWidget(self.label2)
        self.setLayout(layout)

    def initializePage(self):
        path = self.field("version")
        self.label1.setText("raw path is '%s'" % path.toString())
        self.label2.setText("string path is '%s'" % path)

if __name__ == '__main__':
    import sys
    app = QtGui.QApplication(sys.argv)
    wizard = VariantWizard()
    wizard.show()
    sys.exit(app.exec_())