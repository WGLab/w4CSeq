var main = function() {
  $('.mother').click(function() {
    $('.mother').removeClass('current');
    $('.description').hide();

    $(this).addClass('current');
    $(this).children('.description').show();
  });
}

$(document).ready(main);

$(document).ready(function() {
	$('a.mypopover').popover({html:true});
});
function changeTracks()
{
    var insertPos=document.getElementById('dataTrackHere');
    var all=[ {name:'wgEncodeBroadHistoneK562H3k4me1StdPk',cell:'K562',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeSydhTfbsHepg2Tcf7l2UcdPk',cell:'HepG2',experiment:'ChIP-seq_TCF7L2'},
{name:'wgEncodeBroadHistoneHsmmH3k9acStdPk',cell:'HSMM',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHistoneHmecH3k4me1StdPk',cell:'HMEC',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneHepg2H3k27me3StdPk',cell:'HepG2',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeUwDnaseNhlfPkRep1',cell:'NHLF',experiment:'DNase-seq'},
{name:'wgEncodeHaibMethylRrbsHct116StanfordSitesRep1',cell:'HCT116',experiment:'methyl RRBS'},
{name:'tfbsConsSites',cell:'None',experiment:'TFBS Conservation'},
{name:'wgEncodeSydhHistoneHct116H3k27acUcdPk',cell:'HCT116',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneNhlfCtcfStdPk',cell:'NHLF',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeBroadHistoneNhekH3k4me1StdPk',cell:'NHEK',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHmmHepg2HMM',cell:'HepG2',experiment:'chromHMM'},
{name:'wgEncodeBroadHistoneHuvecH3k36me3StdPk',cell:'HUVEC',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeUwTfbsK562CtcfStdPkRep1',cell:'K562',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeOpenChromDnaseHepg2Pk',cell:'HepG2',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneNhlfH3k36me3StdPk',cell:'NHLF',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeBroadHmmK562HMM',cell:'K562',experiment:'chromHMM'},
{name:'wgEncodeUwDnaseHsmmPkRep1',cell:'HSMM',experiment:'DNase-seq'},
{name:'wgEncodeUwHistoneHct116H3k4me3StdPkRep1',cell:'HCT116',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeAwgTfbsUwHelas3CtcfUniPk',cell:'HeLa-S3',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeUwDnaseCaco2PkRep1',cell:'Caco-2',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneH1hescH3k9acStdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHmmHuvecHMM',cell:'HUVEC',experiment:'chromHMM'},
{name:'wgEncodeBroadHistoneH1hescH3k27acStdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneHmecH3k27me3StdPk',cell:'HMEC',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeOpenChromChipHuvecCtcfPk',cell:'HUVEC',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeBroadHistoneHuvecH3k79me2Pk',cell:'HUVEC',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeBroadHmmHsmmHMM',cell:'HSMM',experiment:'chromHMM'},
{name:'wgEncodeUwHistoneCaco2H3k4me3StdPkRep1',cell:'Caco-2',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeAwgTfbsBroadHmecCtcfUniPk',cell:'HMEC',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeBroadHistoneHsmmH3k79me2StdPk',cell:'HSMM',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeSydhHistoneMcf7H3k27acUcdPk',cell:'MCF-7',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeOpenChromFaireHepg2Pk',cell:'HepG2',experiment:'FAIRE-seq'},
{name:'wgEncodeUwDnaseNhekPkRep1',cell:'NHEK',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneH1hescH3k27me3StdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeBroadHistoneHmecH3k9acStdPk',cell:'HMEC',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeOpenChromChipH1hescCtcfPk',cell:'H1-hESC',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeBroadHistoneHsmmH3k27me3StdPk',cell:'HSMM',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeBroadHistoneGm12878H3k9acStdPk',cell:'GM12878',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHistoneHepg2H3k79me2StdPk',cell:'HepG2',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeHaibTfbsHct116CtcfcV0422111PkRep1',cell:'HCT116',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeHaibMethylRrbsH1hescHaibSitesRep1',cell:'H1-hESC',experiment:'methyl RRBS'},
{name:'wgEncodeUwHistoneMcf7H3k4me3StdPkRep1',cell:'MCF-7',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeUwHistoneCaco2H3k36me3StdPkRep1',cell:'Caco-2',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeHaibMethylRrbsA549Dm002p7dHaibSitesRep1',cell:'A549',experiment:'methyl RRBS'},
{name:'wgEncodeBroadHistoneHsmmH3k36me3StdPk',cell:'HSMM',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeHaibMethylRrbsGm12878HaibSitesRep1',cell:'GM12878',experiment:'methyl RRBS'},
{name:'wgEncodeBroadHistoneGm12878H3k36me3StdPk',cell:'GM12878',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeBroadHistoneNhlfH3k27me3StdPk',cell:'NHLF',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeBroadHistoneA549H3k36me3Etoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeBroadHistoneH1hescH3k36me3StdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeUwDnaseHmecPkRep1',cell:'HMEC',experiment:'DNase-seq'},
{name:'wgEncodeHaibMethylRrbsHepg2DukeSitesRep1',cell:'HepG2',experiment:'methyl RRBS'},
{name:'wgEncodeBroadHistoneHsmmH3k9me3StdPk',cell:'HSMM',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeUwTfbsCaco2CtcfStdPkRep1',cell:'Caco-2',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeOpenChromFaireHelas3Pk',cell:'HeLa-S3',experiment:'FAIRE-seq'},
{name:'wgEncodeBroadHistoneA549H3k04me1Etoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHmmNhekHMM',cell:'NHEK',experiment:'chromHMM'},
{name:'wgEncodeBroadHistoneNhekH3k79me2Pk',cell:'NHEK',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeOpenChromDnaseHuvecPk',cell:'HUVEC',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneH1hescH3k4me3StdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeAwgTfbsBroadNhekCtcfUniPk',cell:'NHEK',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeAwgTfbsBroadGm12878CtcfUniPk',cell:'GM12878',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeBroadHistoneHepg2H3k36me3StdPk',cell:'HepG2',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeOpenChromDnaseGm12878Pk',cell:'GM12878',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneHelas3H3k27me3StdPk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeBroadHistoneNhlfH3k4me3StdPk',cell:'NHLF',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeBroadHistoneHepg2H3k27acStdPk',cell:'HepG2',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeUwDnaseHelas3PkRep1',cell:'HeLa-S3',experiment:'DNase-seq'},
{name:'wgEncodeOpenChromFaireH1hescPk',cell:'H1-hESC',experiment:'FAIRE-seq'},
{name:'wgEncodeSydhHistoneHct116H3k04me1UcdPk',cell:'HCT116',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneHuvecH3k4me1StdPk',cell:'HUVEC',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneH1hescH3k4me1StdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeRegTfbsClusteredV3',cell:'None',experiment:'TFBS Region'},
{name:'wgEncodeBroadHistoneK562H3k9me3StdPk',cell:'K562',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeBroadHistoneGm12878H3k27me3StdPk',cell:'GM12878',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeBroadHistoneHuvecH3k4me3StdPk',cell:'HUVEC',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeBroadHistoneGm12878H3k4me1StdPk',cell:'GM12878',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneHelas3H3k79me2StdPk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeBroadHistoneHsmmH3k27acStdPk',cell:'HSMM',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneHuvecH3k09me3Pk',cell:'HUVEC',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeSydhHistoneMcf7H3k36me3bUcdPk',cell:'MCF-7',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeBroadHistoneNhlfH3k79me2Pk',cell:'NHLF',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeHaibMethylRrbsHelas3HaibSitesRep1',cell:'HeLa-S3',experiment:'methyl RRBS'},
{name:'wgEncodeBroadHistoneHepg2H3k04me1StdPk',cell:'HepG2',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneNhekH3k27acStdPk',cell:'NHEK',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneGm12878H3k27acStdPk',cell:'GM12878',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneNhlfH3k27acStdPk',cell:'NHLF',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeAwgTfbsHaibGm12878Tcf3Pcr1xUniPk',cell:'GM12878',experiment:'ChIP-seq_TCF3'},
{name:'wgEncodeUwDnaseHct116PkRep1',cell:'HCT116',experiment:'DNase-seq'},
{name:'wgEncodeSydhHistoneMcf7H3k09me3UcdPk',cell:'MCF-7',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeBroadHistoneA549H3k27acEtoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneHsmmH3k4me1StdPk',cell:'HSMM',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneH1hescH3k09me3StdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeBroadHistoneHmecH3k09me3Pk',cell:'HMEC',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeSydhHistoneMcf7H3k27me3bUcdPk',cell:'MCF-7',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeHaibMethylRrbsK562HaibSitesRep1',cell:'K562',experiment:'methyl RRBS'},
{name:'wgEncodeBroadHistoneA549H3k04me3Etoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeBroadHistoneK562H3k9acStdPk',cell:'K562',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeHaibMethylRrbsCaco2UwSitesRep1',cell:'Caco-2',experiment:'methyl RRBS'},
{name:'wgEncodeSydhTfbsHct116Tcf7l2UcdPk',cell:'HCT116',experiment:'ChIP-seq_TCF7L2'},
{name:'wgEncodeBroadHistoneGm12878H3k9me3StdPk',cell:'GM12878',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeBroadHistoneK562H3k27acStdPk',cell:'K562',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneNhekH3k27me3StdPk',cell:'NHEK',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeAwgTfbsBroadHsmmCtcfUniPk',cell:'HSMM',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeBroadHistoneHuvecH3k27me3StdPk',cell:'HUVEC',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeUwTfbsA549CtcfStdPkRep1',cell:'A549',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeOpenChromDnaseA549Pk',cell:'A549',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneH1hescH3k79me2StdPk',cell:'H1-hESC',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeOpenChromDnaseH1hescPk',cell:'H1-hESC',experiment:'DNase-seq'},
{name:'wgEncodeOpenChromFaireGm12878Pk',cell:'GM12878',experiment:'FAIRE-seq'},
{name:'wgEncodeBroadHistoneK562H3k36me3StdPk',cell:'K562',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeBroadHistoneGm12878H3k79me2StdPk',cell:'GM12878',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeBroadHistoneK562H3k4me3StdPk',cell:'K562',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeUwDnaseMcf7PkRep1',cell:'MCF-7',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneA549H3k79me2Etoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeBroadHistoneHelas3H3k9acStdPk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHistoneA549H3k09me3Etoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeBroadHistoneHelas3H3k27acStdPk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeOpenChromChipHepg2CtcfPk',cell:'HepG2',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeOpenChromFaireHuvecPk',cell:'HUVEC',experiment:'FAIRE-seq'},
{name:'wgEncodeBroadHistoneNhekH3k9acStdPk',cell:'NHEK',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHistoneA549H3k27me3Etoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeOpenChromFaireA549Pk',cell:'A549',experiment:'FAIRE-seq'},
{name:'wgEncodeBroadHistoneHmecH3k4me3StdPk',cell:'HMEC',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeUwHistoneCaco2H3k27me3StdPkRep1',cell:'Caco-2',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeBroadHistoneHuvecH3k27acStdPk',cell:'HUVEC',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneHmecH3k79me2Pk',cell:'HMEC',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeBroadHistoneNhekH3k09me3Pk',cell:'NHEK',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeHaibMethylRrbsMcf7DukeSitesRep1',cell:'MCF-7',experiment:'methyl RRBS'},
{name:'wgEncodeBroadHistoneHepg2H3k9acStdPk',cell:'HepG2',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHmmGm12878HMM',cell:'GM12878',experiment:'chromHMM'},
{name:'wgEncodeBroadHistoneHepg2H3k09me3Pk',cell:'HepG2',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeBroadHmmNhlfHMM',cell:'NHLF',experiment:'chromHMM'},
{name:'wgEncodeBroadHistoneHelas3H3k04me1StdPk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneHelas3H3k36me3StdPk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeBroadHistoneNhekH3k4me3StdPk',cell:'NHEK',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeSydhTfbsHelas3Tcf7l2UcdPk',cell:'HeLa-S3',experiment:'ChIP-seq_TCF7L2'},
{name:'wgEncodeBroadHistoneHmecH3k36me3StdPk',cell:'HMEC',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeBroadHistoneNhlfH3k9acStdPk',cell:'NHLF',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHistoneNhekH3k36me3StdPk',cell:'NHEK',experiment:'ChIP-seq_H3K36me3'},
{name:'wgEncodeUwTfbsMcf7CtcfStdPkRep1',cell:'MCF-7',experiment:'ChIP-seq_CTCF'},
{name:'wgEncodeBroadHistoneGm12878H3k4me3StdPk',cell:'GM12878',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeBroadHmmH1hescHMM',cell:'H1-hESC',experiment:'chromHMM'},
{name:'wgEncodeBroadHistoneK562H3k27me3StdPk',cell:'K562',experiment:'ChIP-seq_H3K27me3'},
{name:'wgEncodeBroadHistoneHelas3H3k4me3StdPk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeBroadHmmHmecHMM',cell:'HMEC',experiment:'chromHMM'},
{name:'wgEncodeBroadHistoneNhlfH3k4me1StdPk',cell:'NHLF',experiment:'ChIP-seq_H3K4me1'},
{name:'wgEncodeBroadHistoneHelas3H3k09me3Pk',cell:'HeLa-S3',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeBroadHistoneHmecH3k27acStdPk',cell:'HMEC',experiment:'ChIP-seq_H3K27ac'},
{name:'wgEncodeBroadHistoneHsmmH3k4me3StdPk',cell:'HSMM',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeBroadHistoneHepg2H3k4me3StdPk',cell:'HepG2',experiment:'ChIP-seq_H3K4me3'},
{name:'wgEncodeBroadHistoneK562H3k79me2StdPk',cell:'K562',experiment:'ChIP-seq_H3K79me2'},
{name:'wgEncodeUwDnaseK562PkRep1',cell:'K562',experiment:'DNase-seq'},
{name:'wgEncodeBroadHistoneA549H3k09acEtoh02Pk',cell:'A549',experiment:'ChIP-seq_H3K9ac'},
{name:'wgEncodeBroadHistoneNhlfH3k09me3Pk',cell:'NHLF',experiment:'ChIP-seq_H3K9me3'},
{name:'wgEncodeOpenChromFaireNhekPk',cell:'NHEK',experiment:'FAIRE-seq'},
{name:'wgEncodeBroadHistoneHuvecH3k9acStdPk',cell:'HUVEC',experiment:'ChIP-seq_H3K9ac'} ];

    var cell=document.getElementsByClassName('cell');
    var experiment=document.getElementsByClassName('experiment');

    var selectedCell=[];
    var selectedExperiment=[];
    var tracks=[];

    //get selected cell
    for (var i=0;i<cell.length;i++)
    {
	if (cell[i].checked)
	{
	    selectedCell.push(cell[i].value);
	}
    }

    //get selected assay
    for (var i=0;i<experiment.length;i++)
    {
	if (experiment[i].checked)
	{
	    selectedExperiment.push(experiment[i].value);
	}
    }

    //get selected tracks
    for (var i=0;i<all.length;i++)
    {
	if (selectedCell.indexOf(all[i].cell) != -1)
	{
	    if (selectedExperiment.indexOf(all[i].experiment) != -1)
	    {
		tracks.push(all[i].name);
	    }
	}
    }

    //remove all child nodes while saving custom uploads info
    while(insertPos.firstChild)
    {
	insertPos.removeChild(insertPos.firstChild);
    }

    //add selected tracks
    for (var i=0;i<tracks.length;i++)
    {
	if (i>19)
	{
	    alert('At most 20 tracks can be selected.');
	    break;
	}
	var newrow=document.createElement('div');//for new row
	var col=document.createElement('div');//for new checkbox container
	var label=document.createElement('label');
	var checkbox=document.createElement('input');
	newrow.class='row';
	col.class='checkbox';
	label.innerHTML=tracks[i];
	checkbox.type='checkbox';
	checkbox.value=tracks[i];
	checkbox.name='generic_table';
	checkbox.checked=true;

	label.appendChild(checkbox);
	col.appendChild(label);
	newrow.appendChild(col);
	insertPos.appendChild(newrow);
    }

    //add custom tracks upload fields
    for (var i=0;i<20-tracks.length;i++)
    {
	var newrow=document.createElement('div');
	var col=document.createElement('div');
	var span=document.createElement('span');
	var label=document.createElement('label');
	var upload=document.createElement('input');
	newrow.class='row';
	col.class='';
	span.title='Drag or use Ctrl, Shift to select multiple files';
	label.innerHTML='';
	upload.type='file';
	upload.name='custom_table';
	upload.multiple='multiple';

	label.appendChild(upload);
	span.appendChild(label);
	col.appendChild(span);
	newrow.appendChild(col);
	insertPos.appendChild(newrow);
    }
}

function response_to_select_region(caller)
{
    var region_spec_container=$(caller).parentsUntil("div.single_region_spec_head").siblings().filter("div.region_detail_area");
    if($(caller).val()=='snp')
    {
	$(region_spec_container).find('div.snp_region_spec').show();
	$(region_spec_container).children().not('.snp_region_spec').hide();
    }
    else if ($(caller).val()=='gene')
    {
	$(region_spec_container).find('div.gene_region_spec').show();
	$(region_spec_container).children().not('.gene_region_spec').hide();
    }
    else if($(caller).val()=='chr')
    {
	$(region_spec_container).find('div.chr_region_spec').show();
	$(region_spec_container).children().not('.chr_region_spec').hide();
    }
}

function response_to_multi_select_region(caller)
{
    var file_region_spec_container=$("#file_region_specification_div_id");
    var manual_region_spec_container=$("#multi_manual_region_specification_div_id");

    if($(caller).val()=='region_file')
    {
	$(file_region_spec_container).show();
	$(manual_region_spec_container).hide();
    }
    else if($(caller).val()=='multi_region')
    {
	$(manual_region_spec_container).show();
	$(file_region_spec_container).hide();
	$(manual_region_spec_container).find("div.single_region_spec_head").each(
		function()
		{
		  var i=$(this).find("div.region_method_area input:checked");
                      i.trigger("click");
		}
	);
    }
}

//display and hide single region selection and multi-region selection alternatively
function toggle_single_multi_region(caller)
{
    var single_container=$("#region_specification_div_id");
    var multi_container=$("#multi_region_specification_div_id");
    if($('#region_multi_single_button_single_id').prop('checked'))
    {
	$(multi_container).hide();
        $(single_container).show();
	$(single_container).find("div.region_method_area input:checked").trigger("click");
    }
    else if ($('#region_multi_single_button_multi_id').prop('checked'))
    {
	$(single_container).hide();
        $(multi_container).show();
	$(multi_container).find("div.multi_region_method_area input:checked").trigger("click");
    }
}

function clear_datatrack_selection()
{
    var cell=document.getElementsByClassName( 'cell');
    for (var i=0;i<cell.length;i++)
    {
	cell[i].checked=false;
    }
    var experiment=document.getElementsByClassName( 'experiment');
    for (var i=0;i<experiment.length;i++)
    {
	experiment[i].checked=false;
    }
}
function loadExampleSetting()
{
    document.getElementById('qformat_whitespace').checked=true;
    document.getElementById('markercol_id').value='dbSNP135';
    document.getElementById('pvalcol_id').value='P';
    document.getElementById('varAnno_id').value='UChicago_eQTL';
    document.getElementById('source_ref_pop_id').value='1000G_March2012,hg19,EUR';	

    $("#region_multi_single_button_single_id").prop('checked',true);
    $("#region_multi_single_button_multi_id").prop('checked',false);
    toggle_single_multi_region();

    //set snp and flanking, check SNP method, trigger SNP click event
    $("div.region_detail_area input[name='snpflank']").val("20");
    $("div.region_detail_area input[name='refsnp']").val("rs2071278");
    $("div.region_detail_area input[name='snpflank']").
    parentsUntil("div.single_region_spec_head").
    siblings().
    filter("div.region_method_area").
    find("input[value='snp']").
    trigger("click");
    $("div.region_detail_area input[name='snpflank']").
    parentsUntil("div.single_region_spec_head").
    siblings().
    filter("div.region_method_area").
    find("input[value='snp']").
    prop('checked',true);

    document.getElementById('generic_toggle_id').checked=true;
    document.getElementById('anno_toggle_id').checked=true;
    document.getElementById('avinput_id').checked=false;
    document.getElementById('ref_hg19').checked=true;

    //set advanced options
    document.getElementById('ld_toggle_id').checked=true;

    //set interaction options
    $('#interaction_cell_type_k562_id').prop('checked',true);
    $('#interaction_cell_type_gm06690_id').prop('checked',false);
    $('#interaction_type_intra_id').prop('checked',true);
    $('#interaction_type_intra_id').trigger('click');
    $('#interaction_type_inter_id').prop('checked',false);
    $('[name="heatmap_toggle"]').prop('checked',true);


    //clear data track selection
    clear_datatrack_selection();
    document.getElementById('HSMM').checked=true;	
    document.getElementById('K562').checked=true;	
    document.getElementById('HepG2').checked=true;	
    document.getElementById('GM12878').checked=true;	
    document.getElementById('None').checked=true;	

    document.getElementById('TFBS Conservation').checked=true;
    document.getElementById('chromHMM').checked=true;
    document.getElementById('TFBS Region').checked=true;
    changeTracks();

    alert('Example settings loaded.');
}
function loadExampleInput()
{

    $('#query_example_label_div_id').show();
    $('#query_hidden_id').val('1');
    $('#query_file_div_id').hide();

    alert('Example input loaded.');
} 
function hideDetail()
{
    document.getElementById('option_detail_id').style.display='none';
}
function showDetail()
{
    document.getElementById('option_detail_id').style.display='block';
}
function check_before_submission()
{
   if($('#query_hidden_id').val() == '0')
   {//we're NOT using example input
    var query_file=$('#query_file_id').val();
    var query_url=$('#query_URL_id').val();

    if (query_file.length>0 && query_url.length>0)
    {
	alert('Input URL cannot be used with input file!');
	//abort submission
	return false;
    }

    if (query_file.length==0 && query_url.length==0)
    {
	alert('No input!');
	return false;
    }
   }

    //check email

    var email=document.getElementById('email_id').value;
    if (email)
    {
    	var email_pattern=/[\w\-\.]+\@[\w\-\.]+\.[\w\-\.]+/;

    	if (! email_pattern.test(email))
    	{
    	    alert('Illegal email address!');
    	    return false;
    	}
    } else
    {
    //when multiple regions are specified, the job will take a while
    //user should provide email
        if($("#region_multi_single_button_multi_id").prop('checked'))
        {
            alert('Email must be provided when multiple regions are to be plotted');
            return false;
        }
    }

    //check genome build
    var ref=document.getElementById('ref_hg19').checked;
    var ref_ld=document.getElementById('source_ref_pop_id').value.split(',');
    ref_ld=ref_ld[1];

    if (ref)
    {
	ref='hg19';
    } else
    {
	ref='hg18';
    }

    if (ref != ref_ld)
    {
	alert('Reference genome does not match source population');
	return false;
    }

    //check marker column
    var marker=document.getElementById('markercol_id').value;

    if (marker.length==0)
    {
	alert('Marker column name is empty!');
	return false;
    }

    //check P value column
    var p_val=document.getElementById('pvalcol_id').value;

    if (p_val.length==0)
    {
	alert('P value column name is empty!');
	return false;
    }

    //Only letters, numbers, dashes, underscores are allowed in column name
    var col_pat=/[^\w\-]/;

    if (col_pat.test(p_val) || col_pat.test(marker))
    {
	alert('Only letters, numbers, dashes, underscores are allowed in column name.');
	return false;
    }

    //check datatracks
    var generic_toggle=document.getElementById('generic_toggle_id').checked;
    var anno_toggle=document.getElementById('anno_toggle_id').checked;
    var datatrack=$('[name="generic_table"]').filter(':checked');
    var custom_table=document.getElementsByName('custom_table');
    var custom_table_count=0;
    $(custom_table).each(function(){
		    custom_table_count+=this.files.length;
		    });

    if ( generic_toggle || anno_toggle )
    {
	if (datatrack.length==0 && custom_table_count==0)
	{
	    alert ("No annotation data tracks selected or uploaded \nwhile generic plot or annotation is enabled.");
	    return false;
	}
    }
    if ( datatrack.length+custom_table_count>20)
    {
	alert('Too many data tracks (max: 20)');
	return false;
    }


    //region specification
    var region_pat=/[ \$\t\r\n\*\|\?\>\<\'\"\,\;\:\[\]\{\}]/;
    var region_method;
    var return_value=true;
    var at_least_one=false;
    
    //when toggled to multi view, the value of toggler is single and vice versa
    if ($('#region_multi_single_button_multi_id').prop('checked'))
    {
    	region_method=$('#multi_region_specification_div_id').find("input[name^='region_method']:checked");//select input elements with name starting with region_method
    }
    else if ($('#region_multi_single_button_single_id').prop('checked'))
    {
    	region_method=$('#region_specification_div_id').find("input[name^='region_method']:checked");//select input elements with name starting with region_method
    }
    
    //user must either upload a region specification file
    //or specify region for each region specification table

	if(! $("input[name='region_file']").val() )
	{
		$(region_method).each(
				function() {
				var region_spec_container=$(this).parentsUntil("div.single_region_spec_head").siblings().filter("div.region_detail_area");
				var refsnp=$(region_spec_container).find("input[name^='refsnp']").val();
				var snpflank=$(region_spec_container).find("input[name^='snpflank']").val();
				var refgene=$(region_spec_container).find("input[name^='refgene']").val();
				var geneflank=$(region_spec_container).find("input[name^='geneflank']").val();
				var generefsnp=$(region_spec_container).find("input[name^='refsnp_for_gene']").val();
				var chr=$(region_spec_container).find("select[name^='chr']").val();
				var start=$(region_spec_container).find("input[name^='start']").val();
				var end=$(region_spec_container).find("input[name^='end']").val();
				var chrrefsnp=$(region_spec_container).find("input[name^='refsnp_for_chr']").val();

					if ($(this).val() == 'snp')
					{	
						if ( (refsnp.length!=0 && snpflank.length==0) || (refsnp.length==0 && snpflank.length!=0) )
						{
							alert('Both reference SNP and flanking region must be supplied.');
							region_spec_container.addClass('has-error');
							return_value=false;
							return false;//stop each iteration
						} else if (region_pat.test(refsnp) || region_pat.test(snpflank))
						{
							alert("Illegal characters found in reference SNP or flanking regin\nPlease remove $ ' \" { } [ ] \\ > < : ; , * | and tab, space, newline.");
							region_spec_container.addClass('has-error');
							return_value=false;
							return false;
						}	
						if( refsnp.length!=0 && snpflank.length!=0)
						{
							at_least_one=true;
						}
					} else if ($(this).val() == 'gene')
					{
						if ( (refgene.length!=0 && geneflank.length==0) || (refgene.length==0 && geneflank.length!=0) )
						{
							alert('Both reference gene and flanking region must be supplied.');
							region_spec_container.addClass('has-error');
							return_value=false;
							return false;
						} else if (region_pat.test(refgene) || region_pat.test(geneflank) || region_pat.test(generefsnp))
						{
							alert("Illegal characters found in reference SNP or flanking regin\nPlease remove $ ' \" { } [ ] \\ > < : ; , * | and tab, space, newline.");
							region_spec_container.addClass('has-error');
							return_value=false;
							return false;
						}	
						if( refgene.length!=0 && geneflank.length!=0)
						{
							at_least_one=true;
						}
					} else if ($(this).val() == 'chr')
					{
						if ( (start.length!=0 && end.length==0) || (start.length==0 && end.length!=0) )
						{
							alert('Both chromosome start and end must be supplied');
							region_spec_container.addClass('has-error');
							return_value=false;
							return false;
						} else if (region_pat.test(chr) || region_pat.test(start) || region_pat.test(end) || region_pat.test(chrrefsnp))
						{
							alert("Illegal characters found in reference SNP or flanking regin\nPlease remove $ ' \" { } [ ] \\ > < : ; , * | and tab, space, newline.");
							region_spec_container.addClass('has-error');
							return_value=false;
							return false;
						}	
						if( start.length!=0 && end.length!=0)
						{
							at_least_one=true;
						}
					}
			}
		);
                $(document).ready( function() {
                    //clear has-error status when this field gets focus
                    $(document).find("div.region_detail_area input").focus( function() {
                    	   $(document).find("div.region_detail_area").removeClass("has-error");
                         });
                 });
		if(!at_least_one && return_value)
		{
			//when return_value is false, an alert has been issued, no need to give more alerts
			alert('No region is specified!');
		}
		return return_value&&at_least_one;
	} else
	{
		return true;
	}
}