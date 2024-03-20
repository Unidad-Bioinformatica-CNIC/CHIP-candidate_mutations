import pymysql
from TargetPanelDBConnector import TargetPanelDBConnector


class TargetPanelMySQLConnector(TargetPanelDBConnector):
    def __init__(self, config):
        try:
            self.connection = pymysql.connect(host=config.get('MYSQL', 'host'),
                                              user=config.get('MYSQL', 'user'),
                                              passwd=config.get('MYSQL', 'passwd'),
                                              database=config.get('MYSQL', 'db'),
                                              port=config.getint('MYSQL', 'port'))
        except Exception as e:
            print(f"ERROR: Can't connect to MySQL DB : {str(e)}")
            exit(1)

    def __del__(self):
        self.connection.close()

    def confirm(self):
        self.connection.commit()

    def add_gene(self, panel_name, symbol, transcript):
        query = f"INSERT INTO target_panel (panel_name, gene, feature) " \
                f"VALUES ('{panel_name}', '{symbol}', '{transcript}')"
        cursor = self.connection.cursor()
        cursor.execute(query)

    def erase(self, panel_name):
        query = f"DELETE FROM target_panel WHERE panel_name = '{panel_name}'"
        cursor = self.connection.cursor()
        cursor.execute(query)
        self.connection.commit()

    def get(self, panel_name):
        query = f"SELECT GENE, FEATURE FROM target_panel WHERE panel_name = '{panel_name}'"
        cursor = self.connection.cursor()
        cursor.execute(query)
        return cursor.fetchall()
