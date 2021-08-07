DEFAULT_DB_NAME = 'trifl'
DEFAULT_REACTIVE_RESIDUE = 'C'
DEFAULT_CIMAGE_PARSER = 'flatfile'
DEFAULT_SILAC_COMBINE_LEVEL = 'protein'
DEFAULT_DTA_FOLDER_NAME = 'dta'
DEFAULT_REPORT_PREFIX = 'report'


pragmas = [
    ('journal_mode', 'wal'),
    ('cache_size', -1000 * 32)
]
