# GET ID counts in subset files

path = 'REFERENCE_2018_mappedCK/'
filelist = ["Cys2018_pantarget_5962.csv","Cys2018_reactive_1431.csv","Lys2018_pantarget_8115.csv","Lys2018_reactive_4425.csv"]

for i in filelist:
    filename = path + i
    df = pd.read_csv(filename)
    idcol = list(df.ID)
    idset = set(idcol)
    print("list lenght: ", len(idcol))
    print("set length: ", len(idset))