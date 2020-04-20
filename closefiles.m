%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% closefiles.m
%%%%%%%%
%%%%%%%% close files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if savdat == 1
  % close files
  fclose(sparf);
  fclose(ssf);
  fclose(sxf);
  fclose(svf);
  fclose(sAf);
  fclose(sPerf);
  fclose(sPAf);
  fclose(spf);
  fclose(sexpf);
  %if sav_fin == 1 % closed after use
  %  fclose(sif);
  %  fclose(sistf);
  %end
  if save_fin == 1
    fclose(sff);
    fclose(sfstf);
  end
  if save_vsim == 1
    fclose(srf);
    fclose(swf);
    fclose(srhof);
    fclose(sTpf);
    fclose(sftotf);
    if save_vforce == 1
      fclose(sflocf);
      fclose(sfintf);
      fclose(sfinthf);
      fclose(sfintvf);
    end
  end
end

