## Notification
<b>These data files are transformed from the test data of [DHS](https://github.com/ZeBraHack0/DHS/tree/main/data) . They are only used for testing the source codes. </b>

<b>The datasets in our experiments are downloaded from [CAIDA-2016](http://www.caida.org/data/passive/passive_2016_dataset.xml) and [CAIDA-2019](http://www.caida.org/data/passive/passive_2019_dataset.xml). If you want to use CAIDA datasets, please register in [CAIDA](http://www.caida.org/home/) first and then apply for the traces.</b>

## Test Data


These data files are binary files in big-endian. Each 13-byte string is a network five-tuple in the format of (srcIP, dstIP, srcPort, dstPort, protocol). For a 13-byte string "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c", the ASCII code should be:

- srcIP = "\x00\x01\x02\x03"
- dstIP = "\x04\x05\x06\x07"
- srcPort = "\x08\x09"
- dstPort = "\x0a\x0b"
- protocol = "\x0c".