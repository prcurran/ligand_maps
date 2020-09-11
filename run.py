from hotspots.hs_pharmacophore import PharmacophoreModel


model = PharmacophoreModel.from_pdb("2vta", "A")

print(model)