from ase.utils import prnt


class Command:
    logfile = None
    args = None
    db = None

    def log(self, *args, **kwargs):
        prnt(file=self.logfile, *args, **kwargs)

    def add_parser(self, subparser):
        pass

    def get_filename(self, name=None, ext=None):
        if name is None:
            if self.args.tag is None:
                filename = 'ase'
            else:
                filename = self.args.tag
        else:
            if '.' in name:
                name = name.rsplit('.', 1)[0]
            if self.args.tag is None:
                filename = name
            else:
                filename = name + '-' + self.args.tag

        if ext:
            filename += '.' + ext

        return filename

    def run(self, atoms, name):
        pass

    def finalize(self):
        pass
