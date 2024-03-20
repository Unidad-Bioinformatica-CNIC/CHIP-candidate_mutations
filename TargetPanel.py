from TargetPanelDBConnector import TargetPanelDBConnector


class TargetPanel:

    def __init__(self, config, db_connector: TargetPanelDBConnector, name=None):

        self.connector = db_connector
        self.name = config.get('GENERAL', 'target_panel') if name is None else name
        self.panel = {gene: feature for gene, feature in self.connector.get(self.name)}

        if self.panel is None:
            print(f"ERROR: Panel {self.name} doesn't exists.")
            exit(1)

    @classmethod
    def create_from_file(cls, config, db_connector: TargetPanelDBConnector, name, file_path):
        def valid_transcript():
            if transcript_raw.startswith('NM_'):
                return transcript_raw.strip().split('.')[0]
            else:
                print(f"WARNING: {transcript_raw} is not a valid RefSeq transcript_id")
                exit(1)

        try:
            # Borra el panel anterior
            db_connector.erase(name)
            # Crea el panel de targets
            with open(file_path) as f:
                for line in f:
                    (symbol, transcript_raw) = line.split()
                    db_connector.add_gene(name, symbol.strip(), valid_transcript())
            db_connector.confirm()
            return cls(config, db_connector, name)
        except FileNotFoundError:
            print(f"ERROR: Target Panel file {file_path} not found")
            exit(1)
        except Exception as e:
            print(f"ERROR: {str(e)}")
            exit(1)

    def delete(self):
        self.connector.erase(self.name)

    def transcripts(self):
        return [feature for feature in self.panel.values()]

    def genes(self):
        return [symbol for symbol in self.panel.keys()]

    def transcript_by_gene(self, gene):
        return self.panel.get(gene, None)
