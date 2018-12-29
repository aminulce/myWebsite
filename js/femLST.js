function splitter(matrix2,elems){//splits a coordinate matrix into x,y,z values and rearranges based on element order
	var numRows = 7; //6+1  the first point is counted twice to close the triangle with six points
//							    .
	                    //     / \
						//    .   .
						//   /     \
						//  .___.___.
	var numRowElems = numel('row',elems);
	elems = math.subtract(elems,1);//COnverting from one based index to zero based index
	var xVals = math.zeros(numRows,numRowElems);
	var yVals = math.zeros(numRows,numRowElems);
	var zVals = math.zeros(numRows,numRowElems);
	var order = math.matrix([[0,3,1,4,2,5,0]]);
	
	for (var el=0; el<numRowElems;el++)
	{
		var matrix = math.zeros(6,2);//since there are 6 points in an element
		for (var kk=0; kk<6; kk++)
		{
			matrix=math.subset(matrix, math.index(kk,0),
			                                 matrix2.subset(math.index(elems.subset(math.index(el,kk)),0)));
			matrix=math.subset(matrix, math.index(kk,1),
			                                 matrix2.subset(math.index(elems.subset(math.index(el,kk)),1)));
		}

		for (var ii=0; ii<numRows;ii++)
		{
			xVals = math.squeeze(math.subset(xVals,math.index(ii,el),matrix.subset(math.index(order.subset(math.index(0,ii)),0))));
			yVals = math.squeeze(math.subset(yVals,math.index(ii,el),matrix.subset(math.index(order.subset(math.index(0,ii)),1))));
			if (numel('column',matrix)==3)
			{
				zVals = math.squeeze(math.subset(zVals,math.index(ii,el),matrix.subset(math.index(ii,2))));
			}
			
		}
	}
	if (numel('column',matrix2)==2)
	{
		return [xVals, yVals];
	}
	else if (numel('column',matrix2)==3)
	{
		return [xVals, yVals, zVals];
	}
	
}

function nodeIndexFromNodeNumberAndDirections(nodes,directions,dimension){
	var numNodes = nodes.length, node,direction;
	var nodeIndices = math.squeeze(math.zeros(numNodes,1)).valueOf();
	console.log(nodes);
	console.log(nodeIndices);
	
	if(dimension==2 && numNodes !=1)
	{
		for(var ii=0;ii<numNodes;ii++)
		{
			node = nodes[ii];
			direction = directions[ii];
			if (direction==1)
			{
				nodeIndices[ii] = node*2-2;
			}
			else if (direction==2)
			{
				console.log('Here');
				nodeIndices[ii] = node*2-1;
				console.log(ii);
			}
			else
				alert('For a 2D analysis second column on nodes only accept 1 and 2');
		}
	return	nodeIndices;	
	}
	else{
		node = nodes[0];
		direction = directions[0];
		if (direction==1)
		{
			nodeIndices = node*2-2;
		}
		else if (direction==2)
		{
			nodeIndices = node*2-1;
		}
		console.log(nodeIndices);
		return	nodeIndices;
	}
}
function numel(rc,matrix){
	if (rc=='row')
	{
		return matrix._size[0];
	}
	else if (rc=='column')
	{
		return matrix._size[1];
	}
}
function lst(nodes,elems,forceMatrix,BC, E, nu){
var t=0.1,CM,cmMult,nuMat;
nuMat = math.matrix([[1,  nu, 0],[nu, 1,  0],[0,  0,  (1-nu)/2]]);
cmMult = E/(1-math.pow(nu,2));
CM= math.matrix(math.multiply(cmMult,nuMat));
var nattau = math.matrix([[2/3, 1/6, 1/6],
						  [1/6, 2/3, 1/6],
						  [1/6, 1/6, 2/3]]);



function reshapeMat(matrix,row,column)
{
	if (row*column!=matrix._size[0])
	{
		var err = 'dimension mismatch while reshaping';
		return err;
	}
	else
	{
		var matrix2 = math.zeros(row,column);
		var count = 0;
		for (var ii=0; ii<row; ii++)
		{
			for (var jj=0; jj<column; jj++)
			{
				matrix2 = math.subset(matrix2,math.index(ii,jj),matrix.subset(math.index(count++)));
			}
		}
		return matrix2;
	}
}
var errMsg = '';
var en;
if(elems.size().length)//building stiffness matrix
{
	numElements = math.size(elems);
	
	var k=math.zeros([numel('row',nodes)*2,numel('row',nodes)*2]);//This k updates after full calculation of an element
	
	var kElement = math.zeros([12,12]);
	for (var ii=0;ii<parseInt(numElements.subset(math.index(0)));ii++)//Loop over elements
	{
		var kIter=math.zeros([numel('row',nodes)*2,numel('row',nodes)*2]);//this k updates at every point
		en = math.squeeze(elems.subset(math.index(ii,[0,1,2,3,4,5])));
		var x = math.squeeze(nodes.subset(math.index(math.subtract(en,1),0)));
		var y = math.squeeze(nodes.subset(math.index(math.subtract(en,1),1)));
		var ke = math.zeros([12,12]);
		var test=1;
		for (var jj=0;jj<3;jj++)//loop over internal goussian points
		{
			var t1 = nattau.subset(math.index(jj,0));
			var t2 = nattau.subset(math.index(jj,1));
			var t3 = nattau.subset(math.index(jj,2));
			
			var dndt1=math.matrix ([4*t1-1,0,0,4*t2,0,4*t3]);
			var dndt2=math.matrix ([0,4*t2-1,0,4*t1,4*t3,0]);
			var dndt3=math.matrix ([0,0,4*t3-1,0,4*t2,4*t1]);
			
			var J = math.matrix([[1,1,1],
			                     [math.multiply(dndt1,x), math.multiply(dndt2,x), math.multiply(dndt3,x)],
								 [math.multiply(dndt1,y), math.multiply(dndt2,y), math.multiply(dndt3,y)]])
								 
								 
			var Jdet = math.det(J)/2;
			var R = math.matrix([[0,0],[1,0],[0,1]]);
			var P = math.multiply(math.inv(J),R);
			
			var dndx=math.add(math.multiply(dndt1,P.subset(math.index(0,0))),
				     math.multiply(dndt2,P.subset(math.index(1,0))),
					 math.multiply(dndt3,P.subset(math.index(2,0))));
					 
			var dndy=math.add(math.multiply(dndt1,P.subset(math.index(0,1))),
				     math.multiply(dndt2,P.subset(math.index(1,1))),
					 math.multiply(dndt3,P.subset(math.index(2,1))));
					 
			var B = math.matrix([[dndx.subset(math.index(0)),
			                      0,
								  dndx.subset(math.index(1)),
								  0,
								  dndx.subset(math.index(2)),
								  0,
								  dndx.subset(math.index(3)),
								  0,
								  dndx.subset(math.index(4)),
								  0,
								  dndx.subset(math.index(5)),
								  0],
			                     [0,
								  dndy.subset(math.index(0)),
								  0,
								  dndy.subset(math.index(1)),
								  0,
								  dndy.subset(math.index(2)),
								  0,
								  dndy.subset(math.index(3)),
								  0,
								  dndy.subset(math.index(4)),
								  0,
								  dndy.subset(math.index(5))],
								 [dndy.subset(math.index(0)),
								  dndx.subset(math.index(0)),
								  dndy.subset(math.index(1)),
								  dndx.subset(math.index(1)),
								  dndy.subset(math.index(2)),
								  dndx.subset(math.index(2)),
								  dndy.subset(math.index(3)),
								  dndx.subset(math.index(3)),
								  dndy.subset(math.index(4)),
								  dndx.subset(math.index(4)),
								  dndy.subset(math.index(5)),
								  dndx.subset(math.index(5)),]]);
								  
			ke = math.add(ke,math.multiply((1/3),Jdet,0.1,math.transpose(B),CM,B));
		}
		var DOF = math.subtract(math.matrix([en.subset(math.index(0))*2-1,
										     en.subset(math.index(0))*2,
										     en.subset(math.index(1))*2-1,
										     en.subset(math.index(1))*2,
										     en.subset(math.index(2))*2-1,
										     en.subset(math.index(2))*2,
										     en.subset(math.index(3))*2-1,
										     en.subset(math.index(3))*2,
										     en.subset(math.index(4))*2-1,
										     en.subset(math.index(4))*2,
										     en.subset(math.index(5))*2-1,
										     en.subset(math.index(5))*2]),1);
		for (jj=0;jj<12;jj++)//putting the element stiffness matrix into global stiffness matrix
		{
			for (var kk=0;kk<12;kk++)
			{
				kIter=math.subset(kIter, math.index(DOF.subset(math.index(jj)), DOF.subset(math.index(kk))), 
							  ke.subset(math.index(jj,kk)));
			}
		}
		k=math.add(k,kIter);
	}
	


for (ii=0; ii<numel('column',BC);ii++)
{
	for (jj=0; jj<numel('row',nodes)*2;jj++)
	{
		k=math.subset(k, math.index(BC.subset(math.index(0,ii)), jj),0);
		k=math.subset(k, math.index(jj,BC.subset(math.index(0,ii))),0);
	}
	k=math.subset(k, math.index(BC.subset(math.index(0,ii)),BC.subset(math.index(0,ii))),1);
}
var dispMatrix = math.multiply(math.inv(k),forceMatrix);
//Get the multiplier as an input through arguments
var multiplier = 100;
var dispMatrix2 = math.multiply(reshapeMat(dispMatrix,numel('row',nodes),2),multiplier);//reshaping the displacemetric to facilitate adding to the node coordinates

var displacedNodes = math.add(nodes,dispMatrix2);
$('#nano').html('test');
const jsonNodes = JSON.stringify(nodes)
const jsonDisplacedNodes = JSON.stringify(displacedNodes)
return[
	jsonNodes,
	jsonDisplacedNodes];
}
else{
	errMsg += 'No element found';
}
};