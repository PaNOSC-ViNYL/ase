
//
// autocomplete
//

function uniq(a)
{
    return a.sort().filter(function(item, pos, ary) 
    {
        return !pos || item != ary[pos - 1];
    })
}

$(function() 
{
    var availableTags = JSON.parse(sessionStorage.getItem("formula"));
    
    $( "#formula-result" ).autocomplete({
        source: availableTags
    });
});

//
// 
//

function ControlFunction(type, element, value)
{
  if(value === 'Reset')
    ns.SetField(element, '');
  else
    ns.SetField(element, value);
    
  document.getElementById("formula-result").focus(); 
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

class inputCtrl
{
    constructor(key, type, index)
    {
        this.m_key = key;
        this.m_type = type; // ['text', 'COMBO', 'INTERVAL', 'CHECK']
        this.m_access = -1;
        this.m_value = [];
    }

    GetQueryString()
    {
        if(this.m_type === "COMBO" || this.m_type === "text")
        {
            return this.m_key + '=' + this.m_value;
        }
        else if(this.m_type === "CHECK")
        {
            return this.m_key;
        }
    }

    TryAssign(tokenID, token)
    {        
        if(this.m_type === "COMBO" || this.m_type === "text")
        {
            // we want a perfect match here
            var checkstr = token.substr(0, this.m_key.length + 1);
            //var indices = ns.getIndicesOf(this.m_key + '=', token);
				
			if(checkstr === this.m_key + '=')
            {
                this.m_access = tokenID;
                this.m_value = token.substr(this.m_key.length + 1, token.length);

                return true;
            }        
        }
        else if(this.m_type === "CHECK")
        {
            if(this.m_key === token)
            {
                this.m_access = tokenID;
                this.m_value = 1;
                return true;
            }
        }

        return false;
    }

    SetOutput()
    {
        if(this.m_access === -1)
            return;

        var element = this.m_key + '-result';

        if(this.m_type === "text")
			document.getElementById(element).value = this.m_value;
        else if(this.m_type === "CHECK")
			document.getElementById(element).checked = this.m_value;
		else
			document.getElementById(element).innerHTML = this.m_value;
    }
};

var ns = (function()
{
    //
	// local variables in namespace
    //

    // search string
	var m_query;

    // tokens
	var m_qlist;
	var m_recognized;

    // controls
    var m_control = [];

	//
	// public functions
    function Init(query)
	{
		m_query = query;

		// set query in forms
		document.getElementById('searchstr').value = query;
		document.getElementById('searchstr2').value = query;

        //
        // setup supported controls
        SetupControls();

        //
        // retrive tokens from search string
        ExtractTokens();

		//
        // init controls as specified in the search string
		SetLoaded();

        //
        // write query
        UpdateSearchQuery();
	}

    function GetControlID(key)
    {
        for(i=0; i<m_control.length; ++i)
        {
            if(m_control[i].m_key == key)
                return i;
        }
        return -1;
    }
    
    function ClearVariable(index)
    {
        m_qlist.splice(m_control[index].m_access, 1);
		m_recognized.splice(m_control[index].m_access, 1);

		for(i=0;i<m_control.length;++i)
		{
			if(m_control[i].m_access > m_control[index].m_access)
				m_control[i].m_access -= 1;
		}

		m_control[index].m_access = -1;
    }

    function CreateVariable(index, key, value)
    {
        m_qlist.push(key + '=' + value);
        m_recognized.push(index);

		m_control[index].m_value = value;
		m_control[index].m_access = m_qlist.length-1;		
    }

	function SetField(key, value)
	{
        var index = GetControlID(key);

        if(index === -1)
            return;

		if(m_control[index].m_access !== -1)
		{
			// the field exists, update its value
			// if the value is empty we remove the variable
			if(value == '')
                ClearVariable(index);
            else
                m_control[index].m_value = value;
		}
		else
		{
			// field has not been set, create it
            CreateVariable(index, key, value);
		}

		// update GUI if it is not a textfield
        m_control[index].SetOutput();

		UpdateSearchQuery();
	}

	//
	// private functions (stays invisible)

    function SetupControls()
    {
        var key = JSON.parse(sessionStorage.getItem("ctrlKeys"));
        var type = JSON.parse(sessionStorage.getItem("ctrlType"));

        for(i=0; i<key.length; i++)
        {
            m_control.push(new inputCtrl(key[i], type[i], i));
        }

        m_control.push(new inputCtrl('formula', 'text', key.length));

        //console.log(m_control);
    }

    function ExtractTokens()
    {
        m_qlist = [];
		m_recognized = [];

		var seperationIndices = getIndicesOf(",", m_query);
		seperationIndices.unshift(-1);
		seperationIndices.push(m_query.length);
			
		for(i=1; i<seperationIndices.length; i++)
		{
			var ss = m_query.substr(seperationIndices[i-1]+1,seperationIndices[i]-seperationIndices[i-1]-1);

            ss = ss.trim();

			if(ss !== '')
			{
				m_qlist.push(ss);
				m_recognized.push(-1);
			}
		}

        //console.log(m_qlist);
    }

    function SetLoaded()
	{
        //
        // For each control find the first appearence in qlist that matches it and try to assign it

		for(j=0; j<m_control.length; j++)
		{
			for(i=0; i<m_qlist.length; i++)
			{
				var indices = getIndicesOf(m_control[j].m_key, m_qlist[i]);                
				
				if(indices.length === 1)
				{
                    if(m_control[j].TryAssign(i, m_qlist[i]) == true)
                    {
						m_recognized[i] = j;
						break;
					}
				}
			}
		}

        //
        // update output controls on the webpage

		for(i=0; i<m_control.length; ++i)
		{
			m_control[i].SetOutput();			
		}
	}

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
				if(m_control[m_recognized[i]].m_value === '')
				{
					continue;
				}
				query += m_control[m_recognized[i]].GetQueryString();
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

    //
    //
    // stuff visible in namespace object
    return {
        Init : Init,
        getIndicesOf : getIndicesOf,
        SetField : SetField
    };
})()