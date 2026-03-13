function spmdsave(fname, wavenumber_factor,wn_k,lat,lat_deg,z,MATRIX_sparse,D,V)
    save(fname, 'wavenumber_factor','wn_k','lat','lat_deg','z','MATRIX_sparse','D','V', '-v7.3')
end
