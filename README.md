# CellBro (wser)


## Interactive Single-Cell Analysis Browser


### Redis dbs

- Session cache: db 0 (server side user sessions)
- flash_cache: db 1 (server side flash message cache)
- message redis: db 2 (messages and status updates from worker to server)
- taskiq redis stream: db 3 (taskiq internal communication)
- taskiq results: db 4 (taskiq task results)
