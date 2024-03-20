import pymysql


def root_dir(config, sample, workflow):
    """ Get `rootdir` directory """
    return f"{config.get('DIR', 'analysis_dir')}/{sample}/{workflow}"


# DB

def connect_db(config):
    try:
        db = pymysql.connect(host=config.get('MYSQL', 'host'),
                             user=config.get('MYSQL', 'user'),
                             passwd=config.get('MYSQL', 'passwd'),
                             database=config.get('MYSQL', 'db'),
                             port=config.getint('MYSQL', 'port'))
        return db

    except Exception as e:
        print(f"ERROR: Can't connect to DB : {str(e)}")
        exit(1)


def exec_query(db, query):
    """ Execute qury on MySQL"""
    cursor = db.cursor()
    cursor.execute(query)
    return cursor.fetchall()

def exec_query_to_dict(db, query):
    """ Execute qury on MySQL"""
    cursor = db.cursor(pymysql.cursors.DictCursor)
    cursor.execute(query)
    return cursor.fetchall()


