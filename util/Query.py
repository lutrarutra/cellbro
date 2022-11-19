

class Query():
    def __init__(self, keys, proposal_keys=None, n_show_keys=20):
        self.keys = keys
        if proposal_keys is None:
            self.proposal_keys = keys[:n_show_keys]
        self.proposal_keys = proposal_keys
        self.n_show_keys = n_show_keys

    def query(self, query):
        self.proposal_keys = []
        for key in self.keys:
            if query.upper() in key.upper():
                self.proposal_keys.append(key)
                if len(self.proposal_keys) == self.n_show_keys:
                    return