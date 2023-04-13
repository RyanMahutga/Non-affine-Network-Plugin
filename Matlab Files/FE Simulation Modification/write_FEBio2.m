function write_FEBio(Input)
%% for FEBio version 3.0

    %load in variables from Input Structure
    base_filename = Input.base_filename;
    filename = Input.filename;
    b=Input.b;
    noElem = Input.noElem;
    bulk = Input.matrixP;
        
    %Initialize matrices
    logfile=[];
    matPart=[];
    mat=[];
    G=[];
    Gx=[];
    Gy=[];
    Gz=[];
    matAxis=[];
    MEe =[];

    %% load base FEB file

    base =regexp(fileread(strcat(pwd,b,base_filename)),'\n','split');
    base=base';

   
    %% Specify Logfile
%     
%     t_s1=find(contains(base,'</plotfile>')); % find start of output
%     t_e1 = find(contains(base,'</Output>')); % find end of output
% 
%     logfile{1,1}=sprintf('		<logfile>');
%     %write F to text file
%     logfile{2,1}=sprintf(strcat('			<element_data file="',Input.outFtensor,'" data="Fxx;Fxy;Fxz;Fyx;Fyy;Fyz;Fzx;Fzy;Fzz;J" name="Fx,Fxy,Fxz,Fyx,Fy,Fyz,Fzx,Fzy,Fz,J">1:',num2str(noElem),'</element_data>'));
%     %write stress to text file
%     logfile{3,1}=sprintf(strcat('			<element_data file="',Input.outStress,'" data="sx;sxy;sxz;sy;syz;sz" name="sx;sxy;sxz;sxy;sy;syz;sxz;syz;sz">1:',num2str(noElem),'</element_data>'));
%     logfile{4,1}=sprintf('		</logfile>');
%     
    %Pressure - will be controlled by load curve 3 **make sure lc3 is in
    %feb file
    pres = find(contains(base,'<pressure'));
    base{pres,1}=sprintf('				<pressure lc="1">%d</pressure>',Input.P);
    
    %% Seperate Material for each element
t_s3=find(contains(base,'<Elements '),1); % start of elements
t_e3 = find(contains(base,'</Elements>'),1); % end of elements

startline = t_s3+1;
    % Partition the elements
    for i =1:noElem
        Ee{1,1}=sprintf('		<Elements type="hex8"  mat="M%d" name="Part%d">',i,i);
        Ee{2,1}=base{startline,1};
        Ee{3,1}=sprintf('		</Elements>');

        matPart = [matPart;Ee];
        startline=startline+1;
    end
    clear i Ee

%     Split up SolidDomain
   t_s5=find(contains(base,'<MeshDomains>'),1);
   t_e5=find(contains(base,'</MeshDomains>'),1);
%    
% %    sl = t_s5+3;
%  matAxis = sprintf('<MeshDomains>');
    for i =1:noElem
        MEe{1,1}=sprintf('		<SolidDomain name="Part%d" mat="M%d"/>',i,i);

        matAxis = [matAxis;MEe];

        clear MEe
    end

    %% Define the Material
    
    t_s4=find(contains(base,'<Material>'),1); % start of material
    t_e4 = find(contains(base,'</Material>'),1); % end of material

    for i=1:noElem
        % HGO fits for each element
        C_input = Input.Matfits.C(i,:);

        % write material section
        %% neo-Hookean
        count =1;
        Ff{count,1}=sprintf('		<material id="%d" name="M%d" type="PeriodicNetwork">',i,i);
        count=count+1;   
        Ff{count,1}=sprintf('			<netNum>%i</netNum>', i);
        count=count+1;
        Ff{count,1}=sprintf('			<netSolve lc="1">%i</netSolve>', 1);
        count=count+1;
        Ff{count,1}=sprintf('			<netSave lc="1">%i</netSave>', 1);
        count=count+1;
        Ff{count,1}=sprintf('			<E>%d</E>',C_input(1));%
        count=count+1;
        Ff{count,1}=sprintf('			<v>%d</v>',Input.matrixP);
        count=count+1;
        Ff{count,1}=sprintf('		</material>');
        %     end

        mat= [mat;Ff];


        clear Ff
    end

    clear i Ff

    t_s6=find(contains(base,'<ElementData'),1);
     t_e6 = find(contains(base,'</ElementData>'),1); 

     dir_elem=[];
     for i=1:noElem
        str1{1} = sprintf('		<ElementData var="mat_axis" elem_set="Part%d">',i);
		str1{2} = sprintf('	<e lid="1">');
        str1{3} = base{t_s6+1+4*(i-1)+1};
        str1{4} = base{t_s6+1+4*(i-1)+2};
		str1{5}= sprintf('	</e>');
        str1{6} = sprintf('		</ElementData>');

        dir_elem = [dir_elem;str1'];
     end

    %% Combine it all and write new FEB

G= vertcat(base(1:t_s4), mat, base(t_e4:t_s3-1), matPart, base(t_e3+1:t_s5),...
    matAxis, base(t_e5:t_s6-1),dir_elem,base(t_e6+1:end));
    fid = fopen(strcat(pwd,b,filename), 'w');
    fprintf(fid, '%s\n', G{:});
    fclose(fid);

%t_s5),matAxis,base(t_e5:
end