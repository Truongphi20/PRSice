bed_file = "tests/data/EUR.QC.bed"
byte_pos = 0
read_size = (483 + 3)// 4

with open(bed_file, "rb") as f:
    f.seek(byte_pos+3)
    result = f.read(read_size)

print("read_size:", read_size)

print(len(result))
print(bin(int.from_bytes(result, byteorder="little")))