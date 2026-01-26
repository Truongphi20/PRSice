bed_file = "tests/data/EUR.QC.bed"
byte_pos = 0
read_size = (483 + 3)// 4

with open(bed_file, "rb") as f:
    f.seek(byte_pos+3)
    result = f.read(read_size)

print("read_size:", read_size)

print(len(result))
# print(bin(int.from_bytes(result, byteorder="little")))

# Read first 24 samples
geno = []
for i in range(25):
    byte_i = i // 4
    pos_i  = i % 4
    # print(byte_i, pos_i, result[byte_i])
    code = (result[byte_i] >> (pos_i*2)) & 0b11
    geno.append(code)

def plink_code_to_vcf_gt(code):
    if code == 0b00: return "1/1"
    if code == 0b10: return "0/1"
    if code == 0b11: return "0/0"
    return "./."
print([plink_code_to_vcf_gt(g) for g in geno])