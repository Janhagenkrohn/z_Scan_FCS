function photonData = readConfoCor3(filepath)

% Author: Jan-Hagen Krohn, MPI for Biochemistry, 2020

% Based on...
%
% Read_PTU_V1.m  
% Original code: Marcus Sackrow, PicoQUant GmbH, December 2013
% read_OTU_V1 by Omri Bar-Elli and Ron Tenne , February 2017
%
% and
%
% Read_ConfoCor3_Raw.m (part of "PAM" PIE Analysis in MATLAB)
% Don C. Lamb lab, see: doi.org/10.1016/j.bpj.2018.02.035

% Reads a .raw file as written by a Zeiss ConfoCor3 and writes it into a
% struct analogous to what is returned by readPTU.m

% INPUT: filepath   : Name/Path of the .raw file to open

% OUTPUT: photonData: Struct with whatever the function could read out from
%                     the file. Many fields are dummies, as Zeiss .raw
%                     files do not contain a whole lot of information.
%                     Nonempty fields are:
%                       Headers: Struct collecting a rogues' gallery of
%                                metadata
%                       ph_sync: Vector of photon arrival times measured in
%                                sync multiples (TTResult_SyncRate!).
%                       ph_dtime:Dummy variable stating a microtime tag for
%                                each photon. Actually just ones,
%                                meaningless.
%                       ph_channel:Channel for each photon. Actually a
%                                repetition of the same number...
%                       TTResult_SyncRate: Number of syncs per second as
%                                time resolution specifier.
%                       MeasDesc_Resolution: Microtime resolution, but
%                                actually dummy value.
%                       File_CreatingTime: Ideally when the file was
%                                written. Not in the header unfortunately,
%                                so what is actually retrieved is the
%                                file's "last modified" tag.
%                       MeasDesc_AcquisitionTime: Acquition time (or
%                                rather, arrival time of the last photon in
%                                measurement) in ms

photonData = struct('ph_sync', [], ...
    'ph_dtime', [], ...
    'ph_channel', [], ...
    'mark_sync', [], ...
    'mark_chan', [], ...
    'mark_dtime', [], ...
    'TTResult_SyncRate', [], ...
    'MeasDesc_Resolution', [], ...
    'File_CreatingTime', [], ...
    'HWSync_Divider', [], ...
    'MeasDesc_AcquisitionTime', []);

fid=fopen(filepath,'r');

%% ASCII file header processing
% read Tag Head
photonData.Headers.Header = fread(fid, 64, '*char')'; % TagHead.Ident
photonData.Headers.Identifier = fread(fid, 4, 'uint32');    % TagHead.Idx
photonData.Headers.Settings = fread(fid, 4, 'uint32');   % TagHead.Typ
dummy = fread(fid, 8, 'uint32'); % Skip 8 4-byte integers?
photonData.TTResult_SyncRate = photonData.Headers.Settings(4);
photonData.MeasDesc_Resolution = 1;% Resolution in picoseconds!

% File creation/modification date
fileInfo = dir(filepath);
photonData.File_CreatingTime = fileInfo.date;

%% Read the T3 mode event records

T3Record = fread(fid, inf, 'uint32');  % all 32 bits:
fclose(fid);
photonData.ph_sync = cumsum(T3Record);
clear T3Record

% Other photon-wise information
photonData.ph_channel = ones(size(photonData.ph_sync)) * str2double(photonData.Headers.Header(end));
photonData.ph_dtime = ones(size(photonData.ph_sync));

photonData.MeasDesc_AcquisitionTime = photonData.ph_sync(end) * 1E3 / photonData.TTResult_SyncRate;
end
