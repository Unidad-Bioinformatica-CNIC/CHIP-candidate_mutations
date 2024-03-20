from abc import ABC, abstractmethod


class TargetPanelDBConnector(ABC):


    @abstractmethod
    def confirm(self):
        pass

    @abstractmethod
    def add_gene(self, panel_name, symbol, transcript):
        pass

    @abstractmethod
    def erase(self):
        pass

    @abstractmethod
    def get(self):
        pass
