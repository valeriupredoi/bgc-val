from shelve import open as shopen
s = shopen('/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/u-ad980/u-ad980_TotalOMZVolume.shelve')
print s.keys()
modeldata = s['modeldata']
readFiles = s['readFiles']
s.close()
modeldata[('Global', 'Surface', 'sum')].pop(1300.4136986301369,None)
readFiles.remove('/group_workspaces/jasmin2/ukesm/BGC_data/u-ad980/u-ad980o_1y_13001201_13011130_ptrc_T.nc')
s = shopen('/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/u-ad980/u-ad980_TotalOMZVolume.shelve')
print s.keys()
s['modeldata'] = modeldata
s['readFiles'] = readFiles
s.close()

