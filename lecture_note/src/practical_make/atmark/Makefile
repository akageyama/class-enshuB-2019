.PHONY: run clean

run: message3.txt
	cat message3.txt

message1.txt:
	@echo Hello, world. I am $(USER) in `hostname` > message1.txt

message2.txt: message1.txt
	@cat message1.txt | tr [a-z] [A-Z] > message2.txt   # To capital letters

message3.txt: message2.txt
	@cat message2.txt > message3.txt   # Repeat it...
	@cat message2.txt >> message3.txt  # ... for two times.

clean:
	@rm -f message1.txt message2.txt message3.txt
	@rm -f message1.txt~ message2.txt~ message3.txt~

