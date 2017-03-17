
//
// 
//

function ControlFunction(type, element, value)
{
  if(value === 'Reset')
    ns.SetField(element, '');
  else
    ns.SetField(element, value);
    
  document.getElementById("name-result").focus(); 
}
        
//
// collapse well storage
//

$(document).ready(function()
{
  // tooltips
  $('[data-toggle="tooltip"]').tooltip();

  $(".collapse").on('shown.bs.collapse', function()
  {
    if (this.id) 
    {
      sessionStorage[this.id] = 'show';
    }
  }).on('hidden.bs.collapse', function()
  {
    if (this.id) 
    {
      sessionStorage.removeItem(this.id);
    }    
  }).each(function()
  {
    if(this.id && sessionStorage[this.id] === 'show')
    {
        $(this).collapse('show');
    }
  });
})

var ns = (function() 
{
	// local variables in namespace
	var m_query;

	var m_qlist;
	var m_recognized;
	
	const m_supported = ['name', 'xc', 'phase'];
	const m_elementToAccess = {'name':0, 'xc':1, 'phase':2};
	var m_access = [-1, -1, -1];
	var m_value = ['', '', ''];

	//
	// public functions
  function Init(query)
	{
		m_query = query;

		// set query in forms
		document.getElementById('searchstr').value = query;
		document.getElementById('searchstr2').value = query;

		m_qlist = [];
		m_recognized = [];


		var seperationIndices = getIndicesOf(",", m_query);
		seperationIndices.unshift(-1);
		seperationIndices.push(m_query.length);
			
		for(i=1; i<seperationIndices.length; i++)
		{
			var ss = m_query.substr(seperationIndices[i-1]+1,seperationIndices[i]-seperationIndices[i-1]-1);

			if(ss !== '')
			{
				m_qlist.push(ss);			
				m_recognized.push(-1);
			}
		}

		SetLoaded();
	}
    
	function SetField(elementID, value)
	{
		var index = m_elementToAccess[elementID];
		var element = '';
		if(m_access[index] !== -1)
		{
			//the field exists

			m_value[index] = value;

			// if the value is empty we clear the variable
			if(value == '')
			{
				m_qlist.splice(m_access[index], 1);
				m_recognized.splice(m_access[index], 1);

				for(i=0;i<m_access.length;++i)
				{
					if(m_access[i] > m_access[index])
						m_access[i] -= 1;
				}

				m_access[index] = -1;
			}

			element = m_supported[index] + '-result';
		}
		else
		{
			// field has not been set, create

			m_qlist.push(elementID + '=' + value);
			m_value[index] = value;
			m_access[index] = m_qlist.length-1;
			m_recognized.push(index);

			element = elementID + '-result';
		}

		// update GUI if it is not a textfield
		if(m_supported[index] !== 'name')
			document.getElementById(element).innerHTML = m_value[index];

		UpdateSearchQuery();
	}

	//
	// private functions (stays invisible)

	function UpdateSearchQuery()
	{
		var query = '';

		for(i=0; i<m_qlist.length; i++)
		{
			if(m_recognized[i] === -1)
			{
				query += m_qlist[i];
			}
			else
			{
				if(m_value[m_recognized[i]] === '')
				{
					continue;
				}
				query += m_supported[m_recognized[i]] + '=' + m_value[m_recognized[i]];
			}

			if(i !== m_qlist.length-1)
			{
				query += ',';
			}
		}

		document.getElementById('searchstr').value = query;
		document.getElementById('searchstr2').value = query;
	}

    function getIndicesOf(searchStr, str, caseSensitive) 
	{
		var searchStrLen = searchStr.length;
		if (searchStrLen == 0) {
			return [];
		}
		var startIndex = 0, index, indices = [];
		if (!caseSensitive) {
			str = str.toLowerCase();
			searchStr = searchStr.toLowerCase();
		}
		while ((index = str.indexOf(searchStr, startIndex)) > -1) {
			indices.push(index);
			startIndex = index + searchStrLen;
		}
		return indices;
	}

	function SetLoaded()
	{
		for(j=0; j<m_supported.length; j++)
		{
			for(i=0; i<m_qlist.length; i++)
			{
				var indices = getIndicesOf(m_supported[j], m_qlist[i]);
				
				if(indices.length === 1)
				{
					if((index = m_qlist[i].indexOf('=')) > -1)
					{
						m_value[j] = m_qlist[i].substr(index+1, m_qlist[i].length);
						m_access[j] = i;

						m_recognized[i] = j;
						break;
					}
				}				
			}
		}

		for(i=0; i<m_access.length; ++i)
		{
			if(m_access[i] !== -1)
			{
				var element = m_supported[i] + '-result';

				if(m_supported[i] === 'name')
					document.getElementById(element).value = m_value[i];
				else
					document.getElementById(element).innerHTML = m_value[i];
			}			
		}
	}

    // stuff visible in namespace object
    return {
        Init : Init,
        SetField : SetField
    };
})()