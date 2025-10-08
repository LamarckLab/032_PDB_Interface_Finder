import os, math, csv
from collections import defaultdict

# ======= 参数定义 =======
INPUT_PDB  = r"C:\Users\Lamarck\Desktop\input.pdb" # 输入文件路径
OUTPUT_CSV = r"C:\Users\Lamarck\Desktop\interaction.csv" # 输出文件路径
CHAIN1 = "A"
CHAIN2 = "B"
CUTOFF = 4.0          # Å
HEAVY_ONLY = True     # 仅计算重原子（忽略氢）
STD_AA_ONLY = True    # 仅考虑标准氨基酸（排除配体/水）

# ======= 标准氨基酸集合 =======
STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL"
}

# ====== 仅处理 “ATOM” 和 “HETATM” 开头的坐标行 ======
def is_atom_line(line: str) -> bool:
    return line.startswith("ATOM") or line.startswith("HETATM")

# ======按照PDB固定列宽进行切片解析  ======
def parse_atom_record(line: str):
    try:
        record = line[0:6].strip()
        atom_name = line[12:16].strip() # 原子名
        altloc = line[16].strip() # 构象标记
        resname = line[17:20].strip()
        chain_id = line[21].strip()
        resseq = line[22:26].strip()
        icode  = line[26].strip()
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        element = line[76:78].strip() if len(line) >= 78 else ""
    except Exception:
        return None
    return {
        "record": record,
        "atom": atom_name,
        "altloc": altloc,
        "resname": resname,
        "chain": chain_id,
        "resseq": resseq,
        "icode": icode,
        "x": x, "y": y, "z": z,
        "element": element
    }

def is_hydrogen(a) -> bool:
    el = (a.get("element") or "").upper()
    if el in ("H", "D"):
        return True
    return a["atom"].upper().startswith("H")

def residue_key(a):
    # 唯一标识一个残基：链、编号、插入码、三字母名
    return (a["chain"], a["resseq"], a["icode"], a["resname"])

def load_chain_atoms(pdb_path, chain_id, heavy_only=True, std_aa_only=True):
    residues = defaultdict(list)  # key -> list[(x,y,z,atom_name)]
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not is_atom_line(line):
                continue
            a = parse_atom_record(line)
            if not a:
                continue
            if a["chain"] != chain_id:
                continue
            # 只保留 altLoc 为空或 'A' 的构象
            if a["altloc"] not in ("", "A"):
                continue
            if heavy_only and is_hydrogen(a):
                continue
            if std_aa_only and a["resname"] not in STANDARD_AA:
                # 排除水、配体、糖等
                continue
            residues[residue_key(a)].append((a["x"], a["y"], a["z"], a["atom"]))
    # 清理空残基
    return {k: v for k, v in residues.items() if v}

def pretty_res_label(key):
    chain, resseq, icode, resname = key
    icode = icode or ""
    return f"{resname}-{resseq}{icode}:{chain}"

def min_distance_between_residues(res_atoms, lig_atoms):
    min_d2 = float("inf")
    for px, py, pz, _ in res_atoms:
        for qx, qy, qz, _ in lig_atoms:
            dx = px - qx; dy = py - qy; dz = pz - qz
            d2 = dx*dx + dy*dy + dz*dz
            if d2 < min_d2:
                min_d2 = d2
    return math.sqrt(min_d2) if min_d2 < float("inf") else float("inf")

def find_contacts(pdb_path, chain1, chain2, cutoff, heavy_only=True, std_aa_only=True):
    chain1_res = load_chain_atoms(pdb_path, chain1, heavy_only, std_aa_only)
    chain2_res = load_chain_atoms(pdb_path, chain2, heavy_only, std_aa_only)

    if not chain1_res:
        print(f"[WARN] Chain {chain1} 没有可用残基（可能被过滤或不存在）")
    if not chain2_res:
        print(f"[WARN] Chain {chain2} 没有可用残基（可能被过滤或不存在）")

    cutoff2 = cutoff * cutoff
    contacts = []  # list[(key1, key2, min_dist)]

    # 朴素 O(N*M) 遍历；对于典型 RBD-ACE2 接口已经足够
    for k1, atoms1 in chain1_res.items():
        for k2, atoms2 in chain2_res.items():
            hit = False
            # 先做一次“存在性”快速筛选：只要有一对原子 < cutoff 即命中
            for px, py, pz, _ in atoms1:
                for qx, qy, qz, _ in atoms2:
                    dx = px - qx; dy = py - qy; dz = pz - qz
                    if dx*dx + dy*dy + dz*dz < cutoff2:
                        hit = True
                        break
                if hit:
                    break
            if hit:
                mind = min_distance_between_residues(atoms1, atoms2)
                contacts.append((k1, k2, mind))

    contacts.sort(key=lambda x: x[2])  # 按最小距离升序
    return contacts

def save_contacts_csv(contacts, out_csv):
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["chain1_resname","chain1_resseq","chain1_icode",
                    "chain2_resname","chain2_resseq","chain2_icode",
                    "min_distance_A"])
        for k1, k2, d in contacts:
            c1, r1, i1, n1 = k1  # chain, resseq, icode, resname
            c2, r2, i2, n2 = k2
            w.writerow([n1, r1, i1, n2, r2, i2, f"{d:.3f}"])

def main():
    if not os.path.exists(INPUT_PDB):
        print(f"[ERROR] 找不到输入文件：{INPUT_PDB}")
        return

    print("== 接触计算开始 ==")
    print(f"PDB路径 : {INPUT_PDB}")
    print(f"链对   : {CHAIN1} vs {CHAIN2}")
    print(f"阈值   : {CUTOFF:.2f} Å")
    print(f"重原子 : {HEAVY_ONLY} | 仅标准氨基酸 : {STD_AA_ONLY}")

    contacts = find_contacts(
        INPUT_PDB, CHAIN1, CHAIN2, CUTOFF,
        heavy_only=HEAVY_ONLY, std_aa_only=STD_AA_ONLY
    )
    save_contacts_csv(contacts, OUTPUT_CSV)

    print(f"\n共找到可能的结合残基对：{len(contacts)}")
    preview_n = min(10, len(contacts))
    if preview_n > 0:
        print("\n示例（前10条或更少）：")
        print("# Residue(chain1)\tResidue(chain2)\tMinDist(Å)")
        for k1, k2, d in contacts[:preview_n]:
            print(f"{pretty_res_label(k1)}\t{pretty_res_label(k2)}\t{d:.3f}")

    print(f"\n已保存 CSV：{OUTPUT_CSV}")
    print("== 完成 ==")

if __name__ == "__main__":
    main()




