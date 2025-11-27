# top level Makefile
# this "simply" runs the make in the two subdirectories: odelib and src


ODELIB =  $(PWD)/SolverLib
export ODELIB

SUBDIRS = src 

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS): 
	make -C $@

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean ;\
	done

