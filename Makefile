.PHONY: debug


OVERRIDE_FILE ?= compose.override.yaml
ifneq ($(wildcard $(OVERRIDE_FILE)),)
    OVERRIDE_FLAG := -f $(OVERRIDE_FILE)
else
    OVERRIDE_FLAG :=
endif

ENV_FILE ?= .env
ifneq ($(wildcard $(ENV_FILE)),)
	ENV_FILE_FLAG := --env-file $(ENV_FILE)
else
	ENV_FILE_FLAG :=
endif


COMPOSE_DEV := docker compose -f compose.yaml $(OVERRIDE_FLAG) -p cellbro-dev $(ENV_FILE_FLAG)

debug:
	@echo "Debugging..."
	@echo "Current directory: $(CURDIR)"
	$(COMPOSE_DEV) up -d --remove-orphans
	$(COMPOSE_DEV) logs -f $(LOGS) --tail=100 cellbro-app
