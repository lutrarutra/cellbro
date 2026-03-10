from typing import Union
import sqlalchemy as sa
from sqlalchemy import orm

from .session import SyncSession

class DBHandler:
    Session: orm.scoped_session
    
    def __init__(
        self,
        expire_on_commit: bool = False,
        auto_commit: bool = False
    ):
        self._session: SyncSession | None = None
        self._connection: sa.engine.Connection | None = None
        self.expire_on_commit = expire_on_commit
        self.auto_commit = auto_commit

    def connect(
        self, user: str, password: str, host: str, db: str = "opengsync_db", port: Union[str, int] = 5432
    ) -> None:
        self._url = f"postgresql+psycopg://{user}:{password}@{host}:{port}/{db}"
        self.public_url = f"{self._url.split(':')[0]}://{host}:{port}/{db}"
        self._engine = sa.create_engine(self._url)
        try:
            self._connection = self._engine.connect()
        except Exception as e:
            raise Exception(f"Could not connect to DB '{self.public_url}':\n{e}")

        self.session_factory = orm.sessionmaker(
            bind=self._engine, expire_on_commit=self.expire_on_commit,
            class_=SyncSession
        )
        DBHandler.Session = orm.scoped_session(self.session_factory)

    @property
    def session(self) -> SyncSession:
        if self._session is None:
            raise Exception("Session is not open.")
        return self._session
    
    def open_session(self, autoflush: bool = False) -> SyncSession:
        if self._session is not None:
            print("Session is already open")
            return self._session
        self._session = DBHandler.Session(autoflush=autoflush)
        return self._session  # type: ignore

    def close_session(self, commit: bool | None = None, rollback: bool = False) -> bool:
        """ returns True if db was modified """
        modified = False
        if self._session is None:
            print("Session is already closed or was never opened.")
            return False
        
        if commit is None:
            commit = self.auto_commit

        if commit and not rollback:
            if self.needs_commit:
                try:
                    self._session.commit()
                except Exception:
                    print("Commit failed: - rolling back transaction.")
                    self._session.rollback()
                    raise
                modified = True
        elif rollback:
            print("Rolling back transaction...")
            self._session.rollback()
        else:
            if not commit and self.needs_commit:
                print("Session was not committed, but changes were made. This may lead to data loss. Use 'db.commit()', if you want changes to be written to the database.")
        self._session = DBHandler.Session.remove()
        return modified
    
    def close_connection(self) -> None:
        if self._connection is not None:
            self._connection = self._connection.close()

    
    def __del__(self):
        if self._session is not None:
            self.close_session()
        self.close_connection()
        self._engine.dispose()

    def __enter__(self):
        self.open_session()
        return self.session

    def __exit__(self, *_):
        self.close_session()

    @property
    def needs_commit(self) -> bool:
        if self._session is None:
            return False
        
        return bool(self._session.dirty) or bool(self._session.new) or bool(self._session.deleted)