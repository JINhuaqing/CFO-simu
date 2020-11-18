from pathlib import Path

rootDir  = Path(".")
fs = rootDir.glob("*.RData")
for f in fs:
    if not ("000" in f.stem):
        targetName = f.stem[:-2] + "_1000_" + f.stem[-1] + ".RData"
        targetF = rootDir/targetName
        print("="*10)
        print(targetF)
        print(f)
        f.rename(targetF)
