classdef Geometry < GeometryInterface
% Describing your geometry
%
%  In TIGRE the geometry is stored in an structure. To see documentation
%  about geometry, run:
%     
%     doc('TIGRE/Geometry') %dead link
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------
%%  Geometry definition
%
%                  Detector plane, behind
%              |-----------------------------|
%              |                             |
%              |                             |
%              |                             |
%  Centered    |                             |
%    at O      A V    +--------+             |
%              |     /        /|             |
%     A Z      |    /        / |*D           |
%     |        |   +--------+  |             |
%     |        |   |        |  |             |
%     |        |   |     *O |  +             |
%     *--->y   |   |        | /              |
%    /         |   |        |/               |
%   V X        |   +--------+        U       |
%              .--------------------->-------|
%
%            *S
%
%
%% Geometry structure:
%           -nVoxel:        3x1 matrix of number of voxels in the image
%           -sVoxel:        3x1 matrix with the total size in mm of the image
%           -dVoxel:        3x1 matrix with the size of each of the voxels in mm
%           -nDetector:     2x1 matrix of number of voxels in the detector plane
%           -sDetector:     2x1 matrix with the total size in mm of the detector
%           -dDetector:     2x1 matrix with the size of each of the pixels in the detector in mm
%           -DSD:           1x1 scalar value. Distance Source Detector, in mm
%           -DSO:           1x1 scalar value. Distance Source Origin.
%           -offOrigin:     3x1 matrix with the offset in mm of the centre of the image from the origin.
%           -offDetector:   2x1 matrix with the offset in mm of the centre
%           of the detector from the x axis

% Shows Geometry diagram
    properties
        % VARIABLE                                   DESCRIPTION                    UNITS
        %-------------------------------------------------------------------------------------
        % Distances
        DSD                                     % Distance Source Detector      (mm)
        DSO                                     % Distance Source Origin        (mm)
        % Detector parameters
        nDetector            					% number of pixels              (px)
        dDetector           					% size of each pixel            (mm)
        sDetector                               % total size of the detector    (mm)
        % Image parameters
        nVoxel                                  % number of voxels              (vx)
        sVoxel                                  % total size of the image       (mm)
        dVoxel                                  % size of each voxel            (mm)
        % Offsets
        offOrigin                               % Offset of image from origin   (mm)              
        offDetector                             % Offset of Detector            (mm)
                                                % These two can be also defined
                                                % per angle
                                                
        % Auxiliary 
        accuracy=0.5;                           % Variable to define accuracy of
                                                % 'interpolated' projection
                                                % It defines the amoutn of
                                                % samples per voxel.
                                                % Recommended <=0.5             (vx/sample)

        % Optional Parameters
        % There is no need to define these unless you actually need them in your
        % reconstruction


        COR=0;                                  % y direction displacement for 
                                                % centre of rotation
                                                % correction                   (mm)
                                                % This can also be defined per
                                                % angle

        rotDetector=[0;0;0];                    % Rotation of the detector, by 
                                                % X,Y and Z axis respectively. (rad)
                                                % This can also be defined per
                                                % angle        

        mode='cone';                            % Or 'parallel'. Geometry type. 
    end
    
    methods
        function geo_struct = checkGeo(obj, angles)

            check_fieldnames(obj)

            check_image_data(obj)

            check_detector_data(obj)

            obj = check_DSO_and_DSD(obj, angles);

            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Now we know that optional fields are properly written or they would have
            % been flagged before.
            % We need to be explicit here again.
            % geofields_optional={'offOrigin','offDetector','rotDetector','COR',...
            %                     'mode','accuracy'};

            if isfield(obj,'offOrigin')
                check_offOrigin(obj)
                if isequal(size(obj.offOrigin),[3 1])
                    obj.offOrigin=repmat(obj.offOrigin,[1, size(angles,2)]);
                end
            else
                obj.offOrigin=zeros(3,size(angles,2));
            end

            if isfield(obj,'offDetector')
                check_offDetector(obj)
                if isequal(size(obj.offDetector),[2 1])
                    obj.offDetector=repmat(obj.offDetector,[1, size(angles,2)]);
                end
            else
                obj.offDetector=zeros(2,size(angles,2));
            end

            if isfield(obj,'rotDetector')
                check_rotDetector(obj)
                if isequal(size(obj.rotDetector),[3 1])
                    obj.rotDetector=repmat(obj.rotDetector,[1, size(angles,2)]);
                end
            else
                obj.rotDetector=zeros(3,size(angles,2));
            end

            if isfield(obj,'COR')
                check_COR(obj)
                if isequal(size(obj.COR),[1 1])
                    obj.COR=repmat(obj.COR,[1, size(angles,2)]);
                end
            else
                obj.COR=zeros(1,size(angles,2));
            end

            if isfield(obj,'mode')
                check_mod(obj)
            else
                obj.mode='cone';
            end

            if isfield(obj,'accuracy')
                check_accuracy(obj)
                if obj.accuracy>1
                    warning('geo.accuracy too big, you will ignore image information resulting in wrong reconstruction.\n Change geo.accuracy to smaller or equal than 1.')
                end
            else
                obj.accuracy=0.5;
            end
            
            % finally must convert properties into a struct before passing
            % to a mex function
            geo_struct = to_struct(obj);
        end
    end
    
    methods (Access = protected)
        function check_fieldnames(obj)
            
            geofields_mandatory = get_mandatory_fields();
            geofields_optional = get_optional_fields();

            fnames=fieldnames(obj);

            check_for_appropriate_fields(fnames, geofields_mandatory, geofields_optional)

            check_for_missing_fields(fnames, geofields_mandatory)

            check_if_double_type(obj, geofields_mandatory);
    
            function geofields_mandatory = get_mandatory_fields()
                geofields_mandatory={
                                     'nVoxel';
                                     'sVoxel';
                                     'dVoxel';
                                     'nDetector';
                                     'sDetector';
                                     'dDetector';
                                     'DSO';
                                     'DSD'}';
            end

            function geofields_optional = get_optional_fields()
                geofields_optional={
                                    'offOrigin';
                                    'offDetector';
                                    'rotDetector';
                                    'COR';
                                    'mode';
                                    'accuracy'}';
            end

            function check_for_appropriate_fields(fnames, geofields_mandatory,geofields_optional)
                allfields=horzcat(geofields_mandatory,geofields_optional);

                % Find if geo has fields we do not recongize
                unknown=~ismember(fnames,allfields);
                % there must be not unknown variables
                % TODO: Do we really want to enforce this? Perhaps just a warning?
                assert(~sum(unknown),'TIGRE:checkGeo:BadGeometry',['The following fields are not known by TIGRE:\n' strjoin(fnames(unknown)),'\nMake sure you have not misspelled any field or introduced unnecesary fields.'])
            end

            function check_for_missing_fields(fnames, geofields_mandatory)
                mandatory=ismember(geofields_mandatory,fnames);
                assert(sum(mandatory)==length(geofields_mandatory),'TIGRE:checkGeo:BadGeometry',['The following fields are missing:\n' strjoin(geofields_mandatory(~mandatory))])
            end

            function check_if_double_type(geo, geofields_mandatory)
                for ii=1:length(geofields_mandatory)
                    assert(isa(geo.(geofields_mandatory{ii}),'double'),'TIGRE:checkGeo:BadGeometry',['Field geo.', geofields_mandatory{ii},' is not double type.'])
                end
            end
    
        end

        function check_image_data(obj)
            assert(isequal(size(obj.nVoxel),[3 1]),'TIGRE:checkGeo:BadGeometry','geo.nVoxel should be 3x1')
            assert(isequal(obj.nVoxel,round(obj.nVoxel)),'TIGRE:checkGeo:BadGeometry','geo.nVoxel should be a natural number.')

            assert(isequal(size(obj.sVoxel),[3 1]),'TIGRE:checkGeo:BadGeometry','geo.sVoxel should be 3x1')

            assert(isequal(size(obj.dVoxel),[3 1]),'TIGRE:checkGeo:BadGeometry','geo.sVoxel should be 3x1')

            assert(sum(abs(obj.dVoxel.*obj.nVoxel-obj.sVoxel))<1e-6, 'TIGRE:checkGeo:BadGeometry', 'nVoxel*dVoxel is not sVoxel, something is wrong in the numbers')
        end

        function check_detector_data(obj)
            assert(isequal(size(obj.nDetector),[2 1]),'TIGRE:checkGeo:BadGeometry','geo.nDetector should be 2x1')
            assert(isequal(obj.nDetector,round(obj.nDetector)),'TIGRE:checkGeo:BadGeometry','geo.nDetector should be a natural number.')

            assert(isequal(size(obj.sDetector),[2 1]),'TIGRE:checkGeo:BadGeometry','geo.sDetector should be 2x1')

            assert(isequal(size(obj.dDetector),[2 1]),'TIGRE:checkGeo:BadGeometry','geo.sDetector should be 2x1')

            assert(sum(abs(obj.dDetector.*obj.nDetector-obj.sDetector))<1e-6, 'TIGRE:checkGeo:BadGeometry', 'nDetector*dDetector is not sDetector, something is wrong in the numbers')
        end

        function obj = check_DSO_and_DSD(obj, angles)
            assert(isequal(size(obj.DSD),[1 1]) | isequal(size(obj.DSD),[1 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.DSD Should be 1x1 or 1xsize(angles,2)')

            assert(isequal(size(obj.DSO),[1 1]) | isequal(size(obj.DSO),[1 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.DSD Should be 1x1 or 1xsize(angles,2)')
            if isequal(size(obj.DSD),[1 1])
                obj.DSD=repmat(obj.DSD,[1, size(angles,2)]);
            end
            if isequal(size(obj.DSO),[1 1])
                obj.DSO=repmat(obj.DSO,[1, size(angles,2)]);
            end

            assert(all(obj.DSD>=obj.DSO), 'TIGRE:checkGeo:BadGeometry','DSD shoudl be bigger or equal to DSO');
        end

        function check_offOrigin(obj)
           assert(isequal(size(obj.offOrigin),[3 1]) | isequal(size(obj.offOrigin),[3 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.offOrigin Should be 3x1 or 3xsize(angles,2)')
           assert(isa(obj.offOrigin,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.offOrigin is not double type.' )
        end

        function check_offDetector(obj)
            assert(isequal(size(obj.offDetector),[2 1]) | isequal(size(obj.offDetector),[2 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.offDetector Should be 2x1 or 2xsize(angles,2)')
            assert(isa(obj.offDetector,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.offDetector is not double type.' )
        end

        function check_rotDetector(obj)
            assert(isequal(size(obj.rotDetector),[3 1]) | isequal(size(obj.rotDetector),[3 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.rotDetector Should be 3x1 or 3xsize(angles,2)')
            assert(isa(obj.rotDetector,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.rotDetector is not double type.' )
        end

        function check_COR(obj)
            assert(isequal(size(obj.COR),[1 1]) | isequal(size(obj.COR),[1 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.COR Should be 1x1 or 1xsize(angles,2)')
            assert(isa(obj.COR,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.COR is not double type.' )
        end

        function check_mod(obj)
           assert(ischar(obj.mode),'TIGRE:checkGeo:BadGeometry','geo.mode shoudl be a character array');
           assert(strcmp(obj.mode,'cone')|strcmp(obj.mode,'parallel'),'TIGRE:checkGeo:BadGeometry','geo.mode shoudl ''cone'' or ''parallel''')
        end

        function check_accuracy(obj)
            assert(isscalar(obj.accuracy),'TIGRE:checkGeo:BadGeometry','geo.accuracy should be a scalar');
            assert(isa(obj.accuracy,'double'),'TIGRE:checkGeo:BadGeometry','geo.accuracy should be double');
        end
        
    end
end
