import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                              QFileDialog, QPushButton, QLabel, QTableView, QTabWidget,
                              QCheckBox, QComboBox, QGroupBox, QFileDialog)
from PyQt5.QtCore import Qt, QAbstractTableModel
import qpaMain

class DataFrameModel(QAbstractTableModel):
    def __init__(self, data):
        super().__init__()
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])
            else:
                return str(self._data.index[section])
        return None

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.current_data = None
        self.results = None

    def initUI(self):
        self.setWindowTitle('qPCR数据分析平台')
        self.setGeometry(300, 300, 1200, 800)

        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # 主布局
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)

        # 左侧控制面板
        control_panel = QGroupBox("控制面板")
        control_layout = QVBoxLayout()

        # 文件上传
        self.btn_upload = QPushButton('上传数据文件')
        self.btn_upload.clicked.connect(self.open_file)
        control_layout.addWidget(self.btn_upload)

        # 参考基因选择
        self.cmb_ref_gene = QComboBox()
        control_layout.addWidget(QLabel('内参基因:'))
        control_layout.addWidget(self.cmb_ref_gene)

        # 参数设置
        self.chk_outliers = QCheckBox('去除异常值')
        self.chk_outliers.setChecked(True)
        control_layout.addWidget(self.chk_outliers)

        # 分析按钮
        self.btn_analyze = QPushButton('开始分析')
        self.btn_analyze.clicked.connect(self.run_analysis)
        control_layout.addWidget(self.btn_analyze)

        # 下载按钮
        self.btn_download = QPushButton('下载结果')
        self.btn_download.clicked.connect(self.save_results)
        control_layout.addWidget(self.btn_download)

        control_panel.setLayout(control_layout)
        main_layout.addWidget(control_panel)

        # 右侧结果显示区域
        self.tabs = QTabWidget()
        self.tab_raw = QTableView()
        self.tab_results = QTableView()

        self.tabs.addTab(self.tab_raw, "原始数据")
        self.tabs.addTab(self.tab_results, "处理结果")

        main_layout.addWidget(self.tabs)

    def open_file(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "打开数据文件", "", 
                                                  "CSV Files (*.csv);;Text Files (*.txt)", 
                                                  options=options)
        if file_path:
            self.current_data = qpaMain.load_data(file_path)
            self.update_ref_gene_choices()
            self.show_data(self.tab_raw, self.current_data)

    def update_ref_gene_choices(self):
        if self.current_data is not None:
            self.cmb_ref_gene.clear()
            self.cmb_ref_gene.addItems(self.current_data['Target'].unique())

    def show_data(self, widget, data):
        model = DataFrameModel(data)
        widget.setModel(model)

    def run_analysis(self):
        if self.current_data is not None:
            ref_gene = self.cmb_ref_gene.currentText()
            remove_outliers = self.chk_outliers.isChecked()
            self.results = self.add_error_handler(lambda: qpaMain.preprocess(self.current_data, ref_gene, remove_outliers))
            self.update_results_display()

    def update_results_display(self):
        if self.results is not None:
            final_data = self.results.get('final_data')
            if 'final_data' not in self.results:
                print('预处理结果缺少final_data字段')
                return
            self.show_data(self.tab_results, final_data)

    def save_results(self):
        if self.results is not None:
            options = QFileDialog.Options()
            file_path, _ = QFileDialog.getSaveFileName(self, "保存结果文件", "", 
                                                  "CSV Files (*.csv);;Text Files (*.txt)", 
                                                  options=options)
            if file_path:
                final_data = self.results.get('final_data')
                if not final_data.empty:
                    try:
                        final_data.to_csv(file_path, index=False)
                    except KeyError:
                        print('预处理结果缺少final_data字段')
                else:
                    print('预处理结果为空')

    def add_error_handler(self, func):
        try:
            return func()
        except Exception as e:
            print(f'预处理异常: {e}')
            return None

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())