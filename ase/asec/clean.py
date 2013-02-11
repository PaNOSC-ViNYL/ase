    def clean_json_file(self, names=None):
        self.read(skipempty=False)

        n = len(self.data)

        if names:
            for name in names:
                del self.data[name]
        else:
            self.data = dict((key, value) for key, value in self.data.items()
                             if value)

        filename = self.get_filename(ext='json')
        try:
            self.lock.acquire()
            write_json(filename, self.data)
        finally:
            self.lock.release()

        n -= len(self.data)
        self.log('Cleaned', n, ['tasks', 'task'][n == 1])


        parser.add_argument('--clean', action='store_true',
                            help='Remove unfinished tasks from json file.')
