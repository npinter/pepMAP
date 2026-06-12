import logging

from flask import Flask

from .routes import bp
from .services.storage import clear_store_dir, init_store, start_scheduler

APP_VERSION = '1.0.1'


def create_app(config_object='config.Config'):
    app = Flask(__name__)
    app.config.from_object(config_object)
    app.config['APP_VERSION'] = APP_VERSION

    log_level_name = 'DEBUG' if app.debug else 'INFO'
    log_level = getattr(logging, log_level_name.upper(), logging.INFO)
    logging.basicConfig(
        level=log_level,
        format='%(levelname)s %(name)s:%(message)s',
        force=True
    )
    logging.getLogger('werkzeug').setLevel(logging.INFO)

    store = init_store(app)

    if app.config.get('CLEAR_STORE_ON_START'):
        clear_store_dir(app.config['PEPMAP_STORE_DIR'])

    if app.config.get('START_SCHEDULER'):
        start_scheduler(app, store)

    app.register_blueprint(bp)
    return app
