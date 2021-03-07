TESTS = 00 01

copy: ${TESTS:%=test%/.timestamp}

test%/.timestamp:
	cd $(shell dirname $@) && mkdir -p BC SC && python gen_inputs.py && bash run_me.sh

clean:
	rm -fr test*/SC
	rm -fr test*/BC
	rm -f test*/test.png test*/run_me.sh

