.DEFAULT_GOAL := all

LIB_DIR = lib
MOD_DIR = mod
SRC_DIR = src

.PHONY: mods libs prog test


all:
	$(MAKE) -C $(MOD_DIR)
	$(MAKE) -C $(LIB_DIR)
	$(MAKE) -C $(SRC_DIR)

test: all
	py.test -vvv -s
clean:
	$(MAKE) -C $(LIB_DIR) clean
	$(MAKE) -C $(MOD_DIR) clean
	$(MAKE) -C $(SRC_DIR) clean
