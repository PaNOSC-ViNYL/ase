


class Collection:
    def __init__(self, collection_name):
        self.collection_name = collection_name
        
    def ____(self, name):
        if self....:
            self.read()
        return ...
        
    def connect(self):
        filename = __file__ + self.name + '.json'
        return ase.db.connect(filename)
        
    def read(self, name):
        connection = self.connect(read_only=True)
        
        
s22 = Collection('s22')
g2 = Collection('g2')
